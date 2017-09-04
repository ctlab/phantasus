/**
 * Created by dzenkova on 11/18/16.
 */
phantasus.KmeansTool = function () {
};
phantasus.KmeansTool.prototype = {
  toString: function () {
    return "k-means";
  },
  gui: function () {
    // z-score, robust z-score, log2, inverse log2
    return [{
      name: "number_of_clusters",
      type: "text"
    }, {
      name: "replace_NA_with",
      type: "bootstrap-select",
      options: [{
        name: "mean",
        value: "mean"
      }, {
        name: "median",
        value: "median"
      }]
    }];
  },
  execute: function (options) {
    var project = options.project;
    var dataset = project.getSortedFilteredDataset();
    var trueIndices = phantasus.Util.getTrueIndices(dataset);

    var columnIndices = [];
    var rowIndices = [];
    if (options.input.use_selected_only) {
      var selectedDataset = project.getSelectedDataset();
      var selectedIndices = phantasus.Util.getTrueIndices(selectedDataset);
      columnIndices = selectedIndices.columns.length > 0 ? selectedIndices.columns : trueIndices.columns;
      rowIndices = selectedIndices.rows.length > 0 ? selectedIndices.rows : trueIndices.rows;
    }
    else {
      columnIndices = trueIndices.columns;
      rowIndices = trueIndices.rows;
    }

    var number = parseInt(options.input.number_of_clusters);
    if (isNaN(number)) {
      throw new Error("Enter the expected number of clusters");
    }
    var replacena = options.input.replace_NA_with;
    var esPromise = dataset.getESSession();
    esPromise.then(function (essession) {
      var args = {
        es: essession,
        k: number,
        replacena: replacena
      };
      if (columnIndices.length > 0) {
        args.columns = columnIndices;
      }
      if (rowIndices.length > 0) {
        args.rows = rowIndices;
      }
      var req = ocpu.call("performKmeans", args, function (session) {
        session.getObject(function (success) {
          var clusters = JSON.parse(success);

          var v = dataset.getRowMetadata().getByName("clusters");
          if (v == null) {
            v = dataset.getRowMetadata().add("clusters");
          }
          for (var i = 0; i < dataset.getRowCount(); i++) {
            v.setValue(i, clusters[i]);
          }
          v.getProperties().set("phantasus.dataType", "string");
          project.trigger("trackChanged", {
            vectors: [v],
            render: ["color"]
          });
        })
      }, false, "::" + dataset.getESVariable());
      req.fail(function () {
        throw new Error("Kmeans call to OpenCPU failed" + req.responseText);
      });

    });
  }
};