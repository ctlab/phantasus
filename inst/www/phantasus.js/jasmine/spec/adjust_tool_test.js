describe('adjust_tool_test', function () {

  it('log_then_inverse_log', function () {
    var heatmap = new phantasus.HeatMap({
      dataset: new phantasus.Dataset({
        array: [[1, 2], [3, 4]],
        rows: 2,
        columns: 2
      })
    });

    var newHeatMap = new phantasus.AdjustDataTool().execute({
      heatMap: heatmap,
      project: heatmap.getProject(),
      input: {
        log_2: true
      }
    });
    expect(newHeatMap.getProject().getFullDataset()).toBeDatasetValues(
      new phantasus.Dataset({
        array: [[0, 1], [1.584963, 2]],
        rows: 2,
        columns: 2
      }), 0.00001);

    newHeatMap = new phantasus.AdjustDataTool().execute({
      heatMap: newHeatMap,
      project: newHeatMap.getProject(),
      input: {
        inverse_log_2: true
      }
    });

    expect(newHeatMap.getProject().getFullDataset()).toBeDatasetValues(
      new phantasus.Dataset({
        array: [[1, 2], [3, 4]],
        rows: 2,
        columns: 2
      }), 0.00001);

  });
  it('log', function () {
    var heatmap = new phantasus.HeatMap({
      dataset: new phantasus.Dataset({
        array: [[1, 2], [3, 4]],
        rows: 2,
        columns: 2
      })
    });

    var newHeatMap = new phantasus.AdjustDataTool().execute({
      heatMap: heatmap,
      project: heatmap.getProject(),
      input: {
        log_2: true
      }
    });
    expect(newHeatMap.getProject().getFullDataset()).toBeDatasetValues(
      new phantasus.Dataset({
        array: [[0, 1], [1.584963, 2]],
        rows: 2,
        columns: 2
      }), 0.00001);

  });

  it('log(1+x)', function () {
    var heatmap = new phantasus.HeatMap({
      dataset: new phantasus.Dataset({
        array: [[1, 2], [3, 4]],
        rows: 2,
        columns: 2
      })
    });

    var newHeatMap = new phantasus.AdjustDataTool().execute({
      heatMap: heatmap,
      project: heatmap.getProject(),
      input: {
        one_plus_log_2: true
      }
    });

    expect(newHeatMap.getProject().getFullDataset()).toBeDatasetValues(
      new phantasus.Dataset({
        array: [[1, 1.584963], [2, 2.3219280]],
        rows: 2,
        columns: 2
      }), 0.00001);
  });

  it('z-score', function () {
    // v = (v-m)/u
    var heatmap = new phantasus.HeatMap({
      dataset: new phantasus.Dataset({
        array: [[1, 2, 10], [3, 4, 15]],
        rows: 2,
        columns: 3
      })
    });

    // var r_sd = [ 4.932883, 6.658328 ];
    // var r_mean = [ 4.333333, 7.333333 ];
    var newHeatMap = new phantasus.AdjustDataTool().execute({
      heatMap: heatmap,
      project: heatmap.getProject(),
      input: {
        'z-score': true
      }
    });
    expect(newHeatMap.getProject().getFullDataset()).toBeDatasetValues(
      new phantasus.Dataset({
        array: [[-0.6757374, -0.4730162, 1.1487535],
          [-0.6508140, -0.5006262, 1.1514402]],
        rows: 2,
        columns: 3
      }), 0.00001);
  });

  it('robust z-score', function () {
    // v = (v-m)/u
    var heatmap = new phantasus.HeatMap({
      dataset: new phantasus.Dataset({
        array: [[1, 2, 10], [3, 4, 15]],
        rows: 2,
        columns: 3
      })
    });

    // var r_mad = [ 1.4826, 1.4826 ];
    // var r_median = [ 2, 4];
    var newHeatMap = new phantasus.AdjustDataTool().execute({
      heatMap: heatmap,
      project: heatmap.getProject(),
      input: {
        'robust_z-score': true
      }
    });
    expect(newHeatMap.getProject().getFullDataset()).toBeDatasetValues(
      new phantasus.Dataset({
        array: [[-0.6744908, 0.0000000, 5.3959261],
          [-0.6744908, 0.0000000, 7.4193984]],
        rows: 2,
        columns: 3
      }), 0.00001);
  });

  describe('sweep', function () {
    var TEST_ROW_NAME = 'test';
    var tool, dataset, heatmap;

    beforeEach(function () {
      tool = new phantasus.AdjustDataTool();
      dataset = new phantasus.Dataset({
        array: [[1, 2], [3, 4]],
        rows: 2,
        columns: 2,
        esSession: new Promise(function () {})
      });
      var vecRow = dataset.getRowMetadata().add(TEST_ROW_NAME); vecRow.setArray([10,20]);
      vecRow.getProperties().set(phantasus.VectorKeys.DATA_TYPE, 'number');

      var vecColumn = dataset.getColumnMetadata().add(TEST_ROW_NAME); vecColumn.setArray([50,100]);
      vecColumn.getProperties().set(phantasus.VectorKeys.DATA_TYPE, 'number');

      heatmap = new phantasus.HeatMap({
        dataset: dataset
      });
    });

    describe('row', function () {
      beforeEach(function () {
        tool.sweepTarget = [{value: 'row'}];
        tool.sweepRowColumnSelect = [{value: TEST_ROW_NAME}];
      });


      it('divide', function () {
        tool.sweepAction = [{value: 'Divide'}];

        var newHeatMap = tool.execute({
          heatMap: heatmap,
          project: heatmap.getProject(),
          input: {}
        });

        expect(newHeatMap.getProject().getFullDataset()).toBeDatasetValues(
          new phantasus.Dataset({
            array: [[1/10, 2/10],
              [3/20, 4/20]],
            rows: 2,
            columns: 2
          }), 0.1);
      });

      it('subtract', function () {
        tool.sweepAction = [{value: 'Subtract'}];

        var newHeatMap = tool.execute({
          heatMap: heatmap,
          project: heatmap.getProject(),
          input: {}
        });

        expect(newHeatMap.getProject().getFullDataset()).toBeDatasetValues(
          new phantasus.Dataset({
            array: [[1 - 10, 2 - 10],
              [3 - 20, 4 - 20]],
            rows: 2,
            columns: 2
          }), 1);
      });
    });

    describe('column', function () {
      beforeEach(function () {
        tool.sweepTarget = [{value: 'column'}];
        tool.sweepRowColumnSelect = [{value: TEST_ROW_NAME}];
      });


      it('divide', function () {
        tool.sweepAction = [{value: 'Divide'}];

        var newHeatMap = tool.execute({
          heatMap: heatmap,
          project: heatmap.getProject(),
          input: {}
        });

        expect(newHeatMap.getProject().getFullDataset()).toBeDatasetValues(
          new phantasus.Dataset({
            array: [[1/50, 2/100],
              [3/50, 4/100]],
            rows: 2,
            columns: 2
          }), 0.01);
      });

      it('subtract', function () {
        tool.sweepAction = [{value: 'Subtract'}];

        var newHeatMap = tool.execute({
          heatMap: heatmap,
          project: heatmap.getProject(),
          input: {}
        });

        expect(newHeatMap.getProject().getFullDataset()).toBeDatasetValues(
          new phantasus.Dataset({
            array: [[1 - 50, 2 - 100],
              [3 - 50, 4 - 100]],
            rows: 2,
            columns: 2
          }), 1);
      });
    });
  })

});
