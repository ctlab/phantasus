phantasus.HClusterTool = function () {
};
phantasus.HClusterTool.PRECOMPUTED_DIST = 'Matrix values (for a precomputed distance matrix)';
phantasus.HClusterTool.PRECOMPUTED_SIM = 'Matrix values (for a precomputed similarity matrix)';
phantasus.HClusterTool.Functions = [phantasus.Euclidean, phantasus.Jaccard,
  new phantasus.OneMinusFunction(phantasus.Cosine),
  new phantasus.OneMinusFunction(phantasus.KendallsCorrelation),
  new phantasus.OneMinusFunction(phantasus.Pearson),
  new phantasus.OneMinusFunction(phantasus.Spearman),
  phantasus.HClusterTool.PRECOMPUTED_DIST,
  phantasus.HClusterTool.PRECOMPUTED_SIM];
phantasus.HClusterTool.Functions.fromString = function (s) {
  for (var i = 0; i < phantasus.HClusterTool.Functions.length; i++) {
    if (phantasus.HClusterTool.Functions[i].toString() === s) {
      return phantasus.HClusterTool.Functions[i];
    }
  }
  throw new Error(s + ' not found');
};

phantasus.HClusterTool.createLinkageMethod = function (linkageString) {
  var linkageMethod;
  if (linkageString === 'Average') {
    linkageMethod = phantasus.AverageLinkage;
  } else if (linkageString === 'Complete') {
    linkageMethod = phantasus.CompleteLinkage;
  } else if (linkageString === 'Single') {
    linkageMethod = phantasus.SingleLinkage;
  } else {
    throw new Error('Unknown linkage method ' + linkageString);
  }
  return linkageMethod;
};

phantasus.HClusterTool.execute = function (dataset, input) {
  // note: in worker here
  var linkageMethod = phantasus.HClusterTool
    .createLinkageMethod(input.linkage_method);
  var f = phantasus.HClusterTool.Functions.fromString(input.metric);
  if (f === phantasus.HClusterTool.PRECOMPUTED_DIST) {
    f = 0;
  } else if (f === phantasus.HClusterTool.PRECOMPUTED_SIM) {
    f = 1;
  }
  var rows = input.cluster == 'Rows' || input.cluster == 'Rows and columns';
  var columns = input.cluster == 'Columns'
    || input.cluster == 'Rows and columns';
  var doCluster = function (d, groupByFields) {
    return (groupByFields && groupByFields.length > 0) ? new phantasus.HClusterGroupBy(
      d, groupByFields, f, linkageMethod)
      : new phantasus.HCluster(phantasus.HCluster
      .computeDistanceMatrix(d, f), linkageMethod);
  };

  var rowsHcl;
  var columnsHcl;

  if (rows) {
    rowsHcl = doCluster(
      input.selectedColumnsToUseForClusteringRows ? new phantasus.SlicedDatasetView(dataset,
        null, input.selectedColumnsToUseForClusteringRows) : dataset,
      input.group_rows_by);
  }
  if (columns) {
    columnsHcl = doCluster(
      phantasus.DatasetUtil
        .transposedView(input.selectedRowsToUseForClusteringColumns ? new phantasus.SlicedDatasetView(
          dataset, input.selectedRowsToUseForClusteringColumns, null)
          : dataset), input.group_columns_by);

  }
  return {
    rowsHcl: rowsHcl,
    columnsHcl: columnsHcl
  };
};
phantasus.HClusterTool.prototype = {
  toString: function () {
    return 'Hierarchical Clustering';
  },
  init: function (project, form) {
    form.setOptions('group_rows_by', phantasus.MetadataUtil
      .getMetadataNames(project.getFullDataset().getRowMetadata()));
    form
      .setOptions('group_columns_by', phantasus.MetadataUtil
        .getMetadataNames(project.getFullDataset()
          .getColumnMetadata()));
    form.setVisible('group_rows_by', false);
    form
      .setVisible('cluster_rows_in_space_of_selected_columns_only',
        false);
    form.$form.find('[name=cluster]').on(
      'change',
      function (e) {
        var val = $(this).val();
        var showGroupColumns = false;
        var showGroupRows = false;
        if (val === 'Columns') {
          showGroupColumns = true;
        } else if (val === 'Rows') {
          showGroupRows = true;
        } else {
          showGroupColumns = true;
          showGroupRows = true;
        }
        form.setVisible('group_columns_by', showGroupColumns);
        form.setVisible('group_rows_by', showGroupRows);
        form.setVisible(
          'cluster_columns_in_space_of_selected_rows_only',
          showGroupColumns);
        form.setVisible(
          'cluster_rows_in_space_of_selected_columns_only',
          showGroupRows);
      });
  },
  gui: function () {
    return [{
      name: 'metric',
      options: phantasus.HClusterTool.Functions,
      value: phantasus.HClusterTool.Functions[4].toString(),
      type: 'select'
    }, {
      name: 'linkage_method',
      options: ['Average', 'Complete', 'Single'],
      value: 'Average',
      type: 'select'
    }, {
      name: 'cluster',
      options: ['Columns', 'Rows', 'Rows and columns'],
      value: 'Columns',
      type: 'select'
    }, {
      name: 'group_columns_by',
      options: [],
      type: 'bootstrap-select',
      multiple: true
    }, {
      name: 'group_rows_by',
      options: [],
      type: 'bootstrap-select',
      multiple: true
    }, {
      name: 'cluster_columns_in_space_of_selected_rows_only',
      type: 'checkbox'
    }, {
      name: 'cluster_rows_in_space_of_selected_columns_only',
      type: 'checkbox'
    }];
  },
  execute: function (options) {
    var project = options.project;
    var heatmap = options.heatMap;
    var selectedRowsToUseForClusteringColumns = options.input.cluster_columns_in_space_of_selected_rows_only ? project
      .getRowSelectionModel().getViewIndices().values()
      : null;
    if (selectedRowsToUseForClusteringColumns != null && selectedRowsToUseForClusteringColumns.length === 0) {
      selectedRowsToUseForClusteringColumns = null;
    }
    var selectedColumnsToUseForClusteringRows = options.input.cluster_rows_in_space_of_selected_columns_only ? project
      .getColumnSelectionModel().getViewIndices().values()
      : null;
    if (selectedColumnsToUseForClusteringRows != null && selectedColumnsToUseForClusteringRows.length === 0) {
      selectedColumnsToUseForClusteringRows = null;
    }
    var rows = options.input.cluster == 'Rows'
      || options.input.cluster == 'Rows and columns';
    var columns = options.input.cluster == 'Columns'
      || options.input.cluster == 'Rows and columns';
    options.input.selectedRowsToUseForClusteringColumns = selectedRowsToUseForClusteringColumns;
    options.input.selectedColumnsToUseForClusteringRows = selectedColumnsToUseForClusteringRows;
    var dataset = project.getSortedFilteredDataset();
    if (options.input.background === undefined) {
      options.input.background = true;
    }
    options.input.background = options.input.background && typeof Worker !== 'undefined';
    var rowModelOrder;
    var columnModelOrder;
    if (rows) {
      rowModelOrder = [];
      for (var i = 0; i < dataset.getRowCount(); i++) {
        rowModelOrder[i] = project.convertViewRowIndexToModel(i);
      }
    }
    if (columns) {
      columnModelOrder = [];
      for (var i = 0; i < dataset.getColumnCount(); i++) {
        columnModelOrder[i] = project.convertViewColumnIndexToModel(i);
      }
    }
    if (options.input.background === false) {
      var result = phantasus.HClusterTool.execute(dataset, options.input);
      if (result.rowsHcl) {
        var modelOrder = [];
        for (var i = 0; i < result.rowsHcl.reorderedIndices.length; i++) {
          modelOrder[i] = rowModelOrder[result.rowsHcl.reorderedIndices[i]];
        }
        heatmap.setDendrogram(result.rowsHcl.tree, false,
          modelOrder);
      }
      if (result.columnsHcl) {
        var modelOrder = [];
        for (var i = 0; i < result.columnsHcl.reorderedIndices.length; i++) {
          modelOrder[i] = columnModelOrder[result.columnsHcl.reorderedIndices[i]];
        }
        heatmap.setDendrogram(result.columnsHcl.tree, true, modelOrder);
      }
    } else {
      var subtitle = ['clustering '];
      if (rows) {
        subtitle.push(dataset.getRowCount() + ' row'
          + phantasus.Util.s(dataset.getRowCount()));
      }
      if (columns) {
        subtitle.push(rows ? ', ' : '');
        subtitle.push(dataset.getColumnCount() + ' column'
          + phantasus.Util.s(dataset.getColumnCount()));
      }

      var blob = new Blob(
        ['self.onmessage = function(e) {'
        + 'importScripts(e.data.scripts);'
        + 'self.postMessage(phantasus.HClusterTool.execute(phantasus.Dataset.fromJSON(e.data.dataset), e.data.input));'
        + '}']);

      var url = window.URL.createObjectURL(blob);
      var worker = new Worker(url);

      worker.postMessage({
        scripts: phantasus.Util.getScriptPath(),
        dataset: phantasus.Dataset.toJSON(dataset, {
          columnFields: options.input.group_columns_by || [],
          rowFields: options.input.group_rows_by || [],
          seriesIndices: [0]
        }),
        input: options.input
      });

      worker.onmessage = function (e) {
        var result = e.data;
        if (result.rowsHcl) {
          var modelOrder = [];
          for (var i = 0; i < result.rowsHcl.reorderedIndices.length; i++) {
            modelOrder[i] = rowModelOrder[result.rowsHcl.reorderedIndices[i]];
          }
          heatmap.setDendrogram(result.rowsHcl.tree, false,
            modelOrder);
        }
        if (result.columnsHcl) {
          var modelOrder = [];
          for (var i = 0; i < result.columnsHcl.reorderedIndices.length; i++) {
            modelOrder[i] = columnModelOrder[result.columnsHcl.reorderedIndices[i]];
          }
          heatmap.setDendrogram(result.columnsHcl.tree, true,
            modelOrder);
        }
        worker.terminate();
        window.URL.revokeObjectURL(url);
      };
      return worker;
    }

  }
};
