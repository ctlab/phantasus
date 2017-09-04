phantasus.NearestNeighbors = function () {
};
phantasus.NearestNeighbors.Functions = [phantasus.Cosine, phantasus.Euclidean,
  phantasus.Jaccard, phantasus.KendallsCorrelation, phantasus.Pearson, phantasus.Spearman,
  phantasus.WeightedMean];
phantasus.NearestNeighbors.Functions.fromString = function (s) {
  for (var i = 0; i < phantasus.NearestNeighbors.Functions.length; i++) {
    if (phantasus.NearestNeighbors.Functions[i].toString() === s) {
      return phantasus.NearestNeighbors.Functions[i];
    }
  }
  throw new Error(s + ' not found');
};
phantasus.NearestNeighbors.prototype = {
  toString: function () {
    return 'Nearest Neighbors';
  },
  init: function (project, form) {
    var $selectedOnly = form.$form.find('[name=use_selected_only]')
      .parent();
    form.$form
      .find('[name=compute_nearest_neighbors_of]')
      .on(
        'change',
        function (e) {
          var val = $(this).val();
          if (val === 'selected rows' || val === 'column annotation') {
            $($selectedOnly.contents()[1])
              .replaceWith(
                document
                  .createTextNode(' Use selected columns only'));
          } else {
            $($selectedOnly.contents()[1])
              .replaceWith(
                document
                  .createTextNode(' Use selected rows only'));
          }
          form.setVisible('annotation', false);
          if (val === 'column annotation' || val === 'row annotation') {
            var metadata = val === 'column annotation' ? project.getFullDataset()
              .getColumnMetadata() : project.getFullDataset()
              .getRowMetadata();
            var names = [];
            // get numeric columns only
            for (var i = 0; i < metadata.getMetadataCount(); i++) {
              var v = metadata.get(i);
              if (phantasus.VectorUtil.getDataType(v) === 'number') {
                names.push(v.getName());
              }
            }
            names.sort(function (a, b) {
              a = a.toLowerCase();
              b = b.toLowerCase();
              return (a < b ? -1 : (a === b ? 0 : 1));
            });
            form
              .setOptions('annotation', names);
            form.setVisible('annotation', true);
          }
        });
    $($selectedOnly.contents()[1]).replaceWith(
      document.createTextNode(' Use selected columns only'));
    form.setVisible('annotation', false);
  },
  gui: function () {
    return [{
      name: 'metric',
      options: phantasus.NearestNeighbors.Functions,
      value: phantasus.Pearson.toString(),
      type: 'select'
    }, {
      name: 'compute_nearest_neighbors_of',
      options: ['selected rows', 'selected columns', 'column annotation', 'row annotation'],
      value: 'selected rows',
      type: 'radio'
    }, {
      name: 'use_selected_only',
      type: 'checkbox'
    }, {
      name: 'annotation',
      type: 'bootstrap-select'
    }];
  },
  execute: function (options) {
    var project = options.project;
    var isColumns = options.input.compute_nearest_neighbors_of == 'selected columns' || options.input.compute_nearest_neighbors_of == 'row annotation';
    var isAnnotation = options.input.compute_nearest_neighbors_of == 'column annotation' || options.input.compute_nearest_neighbors_of == 'row annotation';
    var heatMap = options.heatMap;
    var f = phantasus.NearestNeighbors.Functions
      .fromString(options.input.metric);
    var dataset = project.getSortedFilteredDataset();

    if (isColumns) {
      // compute the nearest neighbors of row, so need to transpose
      dataset = phantasus.DatasetUtil.transposedView(dataset);
    }
    var selectedIndices = (isColumns ? project.getColumnSelectionModel()
      : project.getRowSelectionModel()).getViewIndices().values();
    if (!isAnnotation && selectedIndices.length === 0) {
      throw new Error('No ' + (isColumns ? 'columns' : 'rows')
        + ' selected');
    }
    var spaceIndices = null;
    if (options.input.use_selected_only) {
      spaceIndices = (!isColumns ? project.getColumnSelectionModel()
        : project.getRowSelectionModel()).getViewIndices().values();
      dataset = phantasus.DatasetUtil.slicedView(dataset, null,
        spaceIndices);
    }
    var d1 = phantasus.DatasetUtil
      .slicedView(dataset, selectedIndices, null);
    var list1;
    if (isAnnotation) {
      list1 = dataset.getColumnMetadata().getByName(options.input.annotation);
      if (!list1) {
        throw new Error('No annotation selected.');
      }
    } else {
      if (d1.getRowCount() > 1) {
        // collapse each column in the dataset to a single value
        var columnView = new phantasus.DatasetColumnView(d1);
        var newDataset = new phantasus.Dataset({
          name: '',
          rows: 1,
          columns: d1.getColumnCount()
        });
        for (var j = 0, ncols = d1.getColumnCount(); j < ncols; j++) {
          var v = phantasus.Percentile(columnView.setIndex(j), 50);
          newDataset.setValue(0, j, v);
        }
        d1 = newDataset;
      }
      list1 = new phantasus.DatasetRowView(d1);
    }

    var list2 = new phantasus.DatasetRowView(dataset);
    var values = [];
    var v = dataset.getRowMetadata().getByName(f.toString());
    if (v == null) {
      v = dataset.getRowMetadata().add(f.toString());
    }
    for (var i = 0, size = dataset.getRowCount(); i < size; i++) {
      v.setValue(i, f(list1, list2.setIndex(i)));
    }
    if (!isColumns) {
      project.setRowSortKeys([new phantasus.SortKey(f.toString(),
        phantasus.SortKey.SortOrder.DESCENDING)], true);
    } else {
      project.setColumnSortKeys([new phantasus.SortKey(f.toString(),
        phantasus.SortKey.SortOrder.DESCENDING)], true);
    }
    project.trigger('trackChanged', {
      vectors: [v],
      render: ['text'],
      columns: isColumns
    });
  }
};
