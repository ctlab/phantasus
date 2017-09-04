phantasus.TsneTool = function () {
};

phantasus.TsneTool.execute = function (dataset, input) {
  // note: in worker here
  var matrix = [];
  var rows = input.project == 'Rows';
  if (!rows) {
    dataset = new phantasus.TransposedDatasetView(dataset);
  }
  var N = dataset.getRowCount();
  var f = phantasus.HClusterTool.Functions.fromString(input.metric);
  if (f === phantasus.TsneTool.PRECOMPUTED_DIST) {
    for (var i = 0; i < N; i++) {
      matrix.push([]);
      for (var j = i + 1; j < N; j++) {
        matrix[i][j] = dataset.getValue(i, j);
      }
    }
  } else if (f === phantasus.TsneTool.PRECOMPUTED_SIM) {
    var max = phantasus.DatasetUtil.max(dataset);
    for (var i = 0; i < N; i++) {
      matrix.push([]);
      for (var j = i + 1; j < N; j++) {
        matrix[i][j] = max - dataset.getValue(i, j);
      }
    }
  } else {
    var list1 = new phantasus.DatasetRowView(dataset);
    var list2 = new phantasus.DatasetRowView(dataset);
    for (var i = 0; i < N; i++) {
      matrix.push([]);
      list1.setIndex(i);
      for (var j = i + 1; j < N; j++) {
        var d = f(list1, list2.setIndex(j));
        matrix[i][j] = d;
      }
    }
  }
  var opt = {};
  opt.epsilon = input.epsilon;
  opt.perplexity = input.perplexity;
  opt.dim = 2;
  var tsne = new tsnejs.tSNE(opt);
  tsne.initDataDist(matrix);
  for (var k = 0; k < 1000; k++) {
    tsne.step();
  }
  var Y = tsne.getSolution();
  return {solution: Y};

}
;
phantasus.TsneTool.prototype = {
  toString: function () {
    return 't-SNE';
  },
  init: function (project, form) {

  },
  gui: function () {
    return [{
      name: 'metric',
      options: phantasus.HClusterTool.Functions,
      value: phantasus.HClusterTool.Functions[3].toString(),
      type: 'select'
    }, {
      name: 'project',
      options: ['Columns', 'Rows'],
      value: 'Columns',
      type: 'select'
    }, {
      name: 'epsilon',
      value: '10',
      type: 'text',
      help: 'learning rate'
    }, {
      name: 'perplexity',
      value: '30',
      type: 'text',
      help: 'number of effective nearest neighbors'
    }];
  },
  execute: function (options) {
    var project = options.project;
    var heatMap = options.heatMap;
    var rows = options.input.project == 'Rows';
    var dataset = project.getSortedFilteredDataset();
    options.input.epsilon = parseInt(options.input.epsilon);
    options.input.perplexity = parseInt(options.input.perplexity);
    var blob = new Blob(
      ['self.onmessage = function(e) {'
      + 'e.data.scripts.forEach(function (s) { importScripts(s); });'
      + 'self.postMessage(phantasus.TsneTool.execute(phantasus.Dataset.fromJSON(e.data.dataset), e.data.input));'
      + '}']);

    var url = URL.createObjectURL(blob);
    var worker = new Worker(url);

    worker.postMessage({
      scripts: [phantasus.Util.getScriptPath()],
      dataset: phantasus.Dataset.toJSON(dataset, {
        columnFields: [],
        rowFields: [],
        seriesIndices: [0]
      }),
      input: options.input
    });

    worker.onmessage = function (e) {
      if (rows) {
        dataset = new phantasus.TransposedDatasetView(dataset);
      }
      var result = e.data.solution;

      var newDataset = new phantasus.Dataset({
        name: 't-SNE',
        rows: dataset.getColumnCount(),
        columns: 2
      });

      for (var i = 0; i < result.length; i++) {
        newDataset.setValue(i, 0, result[i][0]);
        newDataset.setValue(i, 1, result[i][1]);
      }
      var idVector = newDataset.getColumnMetadata().add('id');
      idVector.setValue(0, 'P1');
      idVector.setValue(1, 'P2');
      newDataset.setRowMetadata(phantasus.MetadataUtil.shallowCopy(dataset.getColumnMetadata()));
      var min = phantasus.DatasetUtil.min(newDataset);
      var max = phantasus.DatasetUtil.max(newDataset);
      new phantasus.HeatMap({
        inheritFromParentOptions: {transpose: !rows},
        name: 't-SNE',
        dataset: newDataset,
        parent: heatMap,
        columns: [{
          field: 'id',
          display: 'text'
        }],
        colorScheme: {
          type: 'fixed',
          map: [{
            value: min,
            color: colorbrewer.Greens[3][0]
          }, {
            value: max,
            color: colorbrewer.Greens[3][2]
          }]
        }
      });
      worker.terminate();
      window.URL.revokeObjectURL(url);
    };
    return worker;
  }
};
