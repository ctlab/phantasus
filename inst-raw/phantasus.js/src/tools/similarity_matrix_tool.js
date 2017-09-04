phantasus.SimilarityMatrixTool = function () {
};

phantasus.SimilarityMatrixTool.Functions = [phantasus.Euclidean,
  phantasus.Jaccard, phantasus.Cosine, phantasus.KendallsCorrelation, phantasus.Pearson, phantasus.Spearman];
phantasus.SimilarityMatrixTool.Functions.fromString = function (s) {
  for (var i = 0; i < phantasus.SimilarityMatrixTool.Functions.length; i++) {
    if (phantasus.SimilarityMatrixTool.Functions[i].toString() === s) {
      return phantasus.SimilarityMatrixTool.Functions[i];
    }
  }
  throw new Error(s + ' not found');
};
phantasus.SimilarityMatrixTool.execute = function (dataset, input) {
  var isColumnMatrix = input.compute_matrix_for == 'Columns';
  var f = phantasus.SimilarityMatrixTool.Functions.fromString(input.metric);
  return phantasus.HCluster.computeDistanceMatrix(
    isColumnMatrix ? new phantasus.TransposedDatasetView(dataset)
      : dataset, f);
};
phantasus.SimilarityMatrixTool.prototype = {
  toString: function () {
    return 'Similarity Matrix';
  },
  init: function (project, form) {

  },
  gui: function () {
    return [{
      name: 'metric',
      options: phantasus.SimilarityMatrixTool.Functions,
      value: phantasus.SimilarityMatrixTool.Functions[4].toString(),
      type: 'select'
    }, {
      name: 'compute_matrix_for',
      options: ['Columns', 'Rows'],
      value: 'Columns',
      type: 'radio'
    }];
  },
  execute: function (options) {
    var project = options.project;
    var heatMap = options.heatMap;
    var isColumnMatrix = options.input.compute_matrix_for == 'Columns';
    var f = phantasus.SimilarityMatrixTool.Functions
      .fromString(options.input.metric);
    var dataset = project.getSortedFilteredDataset();
    var blob = new Blob(
      ['self.onmessage = function(e) {'
      + 'importScripts(e.data.scripts);'
      + 'self.postMessage(phantasus.SimilarityMatrixTool.execute(phantasus.Dataset.fromJSON(e.data.dataset), e.data.input));'
      + '}']);

    var url = window.URL.createObjectURL(blob);
    var worker = new Worker(url);

    worker.postMessage({
      scripts: phantasus.Util.getScriptPath(),
      dataset: phantasus.Dataset.toJSON(dataset, {
        columnFields: [],
        rowFields: [],
        seriesIndices: [0]
      }),
      input: options.input
    });

    worker.onmessage = function (e) {
      var name = heatMap.getName();
      var matrix = e.data;
      var n = isColumnMatrix ? dataset.getColumnCount() : dataset
        .getRowCount();
      var d = new phantasus.Dataset({
        name: name,
        rows: n,
        columns: n
      });
      // set the diagonal
      var isDistance = f.toString() === phantasus.Euclidean.toString()
        || f.toString() === phantasus.Jaccard.toString();
      for (var i = 1; i < n; i++) {
        for (var j = 0; j < i; j++) {
          var value = matrix[i][j];
          d.setValue(i, j, value);
          d.setValue(j, i, value);
        }
      }
      // no need to set diagonal if not distance as array already
      // initialized to 0
      if (!isDistance) {
        for (var i = 0; i < n; i++) {
          d.setValue(i, i, 1);
        }
      }
      var metadata = isColumnMatrix ? dataset.getColumnMetadata()
        : dataset.getRowMetadata();
      d.rowMetadataModel = phantasus.MetadataUtil.shallowCopy(metadata);
      d.columnMetadataModel = phantasus.MetadataUtil.shallowCopy(metadata);
      var colorScheme;
      if (!isDistance) {
        colorScheme = {
          type: 'fixed',
          map: [{
            value: -1,
            color: 'blue'
          }, {
            value: 0,
            color: 'white'
          }, {
            value: 1,
            color: 'red'
          }]
        };
      } else {
        colorScheme = {
          type: 'fixed',
          map: [{
            value: 0,
            color: 'white'
          }, {
            value: phantasus.DatasetUtil.max(d),
            color: 'red'
          }]
        };
      }
      new phantasus.HeatMap({
        colorScheme: colorScheme,
        name: name,
        dataset: d,
        parent: heatMap,
        inheritFromParentOptions: {
          rows: !isColumnMatrix,
          columns: isColumnMatrix
        }
      });
      worker.terminate();
      window.URL.revokeObjectURL(url);
    };
    return worker;
  }
};
