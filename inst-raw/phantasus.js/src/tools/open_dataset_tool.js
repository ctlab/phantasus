phantasus.OpenDatasetTool = function () {
  this.customUrls = [];
};

phantasus.OpenDatasetTool.fileExtensionPrompt = function (file, callback) {
  var ext = phantasus.Util.getExtension(phantasus.Util.getFileName(file));
  var deferred;
  if (ext === 'seg' || ext === 'segtab') {
    this._promptSegtab(function (regions) {
      callback(regions);
    });

  } else {
    callback(null);
  }

};
phantasus.OpenDatasetTool._promptMaf = function (promptCallback) {
  var formBuilder = new phantasus.FormBuilder();
  formBuilder
    .append({
      name: 'MAF_gene_symbols',
      value: '',
      type: 'textarea',
      required: true,
      help: 'Enter one gene symbol per line to filter genes. Leave blank to show all genes.'
    });
  phantasus.FormBuilder
    .showInModal({
      title: 'Gene Symbols',
      html: formBuilder.$form,
      close: 'OK',
      onClose: function () {
        var text = formBuilder.getValue('MAF_gene_symbols');
        var lines = phantasus.Util.splitOnNewLine(text);
        var mafGeneFilter = new phantasus.Map();
        for (var i = 0, nlines = lines.length, counter = 0; i < nlines; i++) {
          var line = lines[i];
          if (line !== '') {
            mafGeneFilter.set(line, counter++);
          }
        }
        var readOptions = mafGeneFilter.size() > 0 ? {
          mafGeneFilter: mafGeneFilter
        } : null;
        promptCallback(readOptions);
      }
    });
};
phantasus.OpenDatasetTool._promptSegtab = function (promptCallback) {
  var formBuilder = new phantasus.FormBuilder();
  formBuilder
    .append({
      name: 'regions',
      value: '',
      type: 'textarea',
      required: true,
      help: 'Define the regions over which you want to define the CNAs. Enter one region per line. Each line should contain region_id, chromosome, start, and end separated by a tab. Leave blank to use all unique segments in the segtab file as regions.'
    });
  phantasus.FormBuilder
    .showInModal({
      title: 'Regions',
      html: formBuilder.$form,
      close: 'OK',
      onClose: function () {
        var text = formBuilder.getValue('regions');
        var lines = phantasus.Util.splitOnNewLine(text);
        var regions = [];
        var tab = /\t/;
        for (var i = 0, nlines = lines.length, counter = 0; i < nlines; i++) {
          var line = lines[i];

          if (line !== '') {
            var tokens = line.split(tab);
            if (tokens.length >= 4) {
              regions.push({
                id: tokens[0],
                chromosome: tokens[1],
                start: parseInt(tokens[2]),
                end: parseInt(tokens[3])
              });
            }
          }
        }
        var readOptions = regions.length > 0 ? {
          regions: regions
        } : null;
        promptCallback(readOptions);
      }
    });
};
phantasus.OpenDatasetTool.prototype = {
  toString: function () {
    return 'Open Dataset';
  },
  _read: function (options, deferred) {
    var _this = this;
    var project = options.project;
    var heatMap = options.heatMap;
    var file = options.input.file;
    var action = options.input.open_file_action;
    var dataset = project.getSortedFilteredDataset();
    deferred.fail(function (err) {
      var message = ['Error opening ' + phantasus.Util.getFileName(file)
      + '.'];
      if (err.message) {
        message.push('<br />Cause: ');
        message.push(err.message);
      }
      phantasus.FormBuilder.showInModal({
        title: 'Error',
        html: message.join('')
      });
    });
    deferred
      .done(function (newDataset) {

        var extension = phantasus.Util.getExtension(phantasus.Util
          .getFileName(file));
        var filename = phantasus.Util.getBaseFileName(phantasus.Util
          .getFileName(file));
        if (action === 'append' || action === 'append columns') {

          // "append": append rows to current dataset
          var appendRows = action === 'append';
          // rename fields?
          _.each(heatMap.options.rows, function (item) {
            if (item.renameTo) {
              var v = newDataset.getRowMetadata().getByName(
                item.field);
              if (v) {
                v.setName(item.renameTo);
              }
            }
          });
          _.each(heatMap.options.columns, function (item) {
            if (item.renameTo) {
              var v = newDataset.getColumnMetadata()
                .getByName(item.field);
              if (v) {
                v.setName(item.renameTo);
              }
            }
          });

          if (heatMap.options.datasetReady) {
            heatMap.options.datasetReady(newDataset);
          }
          var currentDatasetMetadataNames = phantasus.MetadataUtil
            .getMetadataNames(!appendRows ? dataset
              .getRowMetadata() : dataset
              .getColumnMetadata());
          var newDatasetMetadataNames = phantasus.MetadataUtil
            .getMetadataNames(!appendRows ? newDataset
              .getRowMetadata() : newDataset
              .getColumnMetadata());

          if (currentDatasetMetadataNames.length > 1
            || newDatasetMetadataNames.length > 1) {

            _this
              ._matchAppend(
                newDatasetMetadataNames,
                currentDatasetMetadataNames,
                heatMap,
                function (appendOptions) {
                  heatMap
                    .getProject()
                    .setFullDataset(
                      appendRows ? new phantasus.JoinedDataset(
                        dataset,
                        newDataset,
                        appendOptions.current_dataset_annotation_name,
                        appendOptions.new_dataset_annotation_name)
                        : new phantasus.TransposedDatasetView(
                        new phantasus.JoinedDataset(
                          new phantasus.TransposedDatasetView(
                            dataset),
                          new phantasus.TransposedDatasetView(
                            newDataset),
                          appendOptions.current_dataset_annotation_name,
                          appendOptions.new_dataset_annotation_name)),
                      true);

                  if (heatMap.options.renderReady) {
                    heatMap.options
                      .renderReady(heatMap);
                    heatMap.updateDataset();
                  }
                  if (appendRows) {
                    heatMap
                      .getHeatMapElementComponent()
                      .getColorScheme()
                      .setSeparateColorSchemeForRowMetadataField(
                        'Source');

                    var sourcesSet = phantasus.VectorUtil
                      .getSet(heatMap
                        .getProject()
                        .getFullDataset()
                        .getRowMetadata()
                        .getByName(
                          'Source'));
                    sourcesSet
                      .forEach(function (source) {
                        heatMap
                          .autoDisplay({
                            extension: phantasus.Util
                              .getExtension(source),
                            filename: source
                          });
                      });
                  }

                  heatMap.tabManager
                    .setTabTitle(
                      heatMap.tabId,
                      heatMap
                        .getProject()
                        .getFullDataset()
                        .getRowCount()
                      + ' row'
                      + phantasus.Util
                        .s(heatMap
                          .getProject()
                          .getFullDataset()
                          .getRowCount())
                      + ' x '
                      + heatMap
                        .getProject()
                        .getFullDataset()
                        .getColumnCount()
                      + ' column'
                      + phantasus.Util
                        .s(heatMap
                          .getProject()
                          .getFullDataset()
                          .getColumnCount()));
                  heatMap.revalidate();
                });
          } else { // no need to prompt
            heatMap
              .getProject()
              .setFullDataset(
                appendRows ? new phantasus.JoinedDataset(
                  dataset,
                  newDataset,
                  currentDatasetMetadataNames[0],
                  newDatasetMetadataNames[0])
                  : new phantasus.TransposedDatasetView(
                  new phantasus.JoinedDataset(
                    new phantasus.TransposedDatasetView(
                      dataset),
                    new phantasus.TransposedDatasetView(
                      newDataset),
                    currentDatasetMetadataNames[0],
                    newDatasetMetadataNames[0])),
                true);
            if (heatMap.options.renderReady) {
              heatMap.options.renderReady(heatMap);
              heatMap.updateDataset();
            }
            if (appendRows) {
              heatMap
                .getHeatMapElementComponent()
                .getColorScheme()
                .setSeparateColorSchemeForRowMetadataField(
                  'Source');
              var sourcesSet = phantasus.VectorUtil
                .getSet(heatMap.getProject()
                  .getFullDataset()
                  .getRowMetadata().getByName(
                    'Source'));
              sourcesSet.forEach(function (source) {
                heatMap.autoDisplay({
                  extension: phantasus.Util
                    .getExtension(source),
                  filename: source
                });
              });
            }
            heatMap.tabManager.setTabTitle(heatMap.tabId,
              heatMap.getProject().getFullDataset()
                .getRowCount()
              + ' row'
              + phantasus.Util.s(heatMap
                .getProject()
                .getFullDataset()
                .getRowCount())
              + ' x '
              + heatMap.getProject()
                .getFullDataset()
                .getColumnCount()
              + ' column'
              + phantasus.Util.s(heatMap
                .getProject()
                .getFullDataset()
                .getColumnCount()));
            heatMap.revalidate();
          }

        } else if (action === 'overlay') {
          _this
            ._matchOverlay(
              phantasus.MetadataUtil
                .getMetadataNames(newDataset
                  .getColumnMetadata()),
              phantasus.MetadataUtil
                .getMetadataNames(dataset
                  .getColumnMetadata()),
              phantasus.MetadataUtil
                .getMetadataNames(newDataset
                  .getRowMetadata()),
              phantasus.MetadataUtil
                .getMetadataNames(dataset
                  .getRowMetadata()),
              heatMap,
              function (appendOptions) {
                phantasus.DatasetUtil.overlay({
                  dataset: dataset,
                  newDataset: newDataset,
                  rowAnnotationName: appendOptions.current_dataset_row_annotation_name,
                  newRowAnnotationName: appendOptions.new_dataset_row_annotation_name,
                  columnAnnotationName: appendOptions.current_dataset_column_annotation_name,
                  newColumnAnnotationName: appendOptions.new_dataset_column_annotation_name
                });
              });
        } else if (action === 'open') { // new tab
          if (newDataset.length && newDataset.length > 0) {
            for (var i = 0; i < newDataset.length; i++) {
              new phantasus.HeatMap({
                dataset: newDataset[i],
                name: newDataset[i].seriesNames[0],
                parent: heatMap,
                inheritFromParent: false
              });
            }
          } else {
            new phantasus.HeatMap({
              dataset: newDataset,
              parent: heatMap,
              inheritFromParent: false
            });
          }

        } else {
          // console.log('Unknown action: ' + action);
        }

        heatMap.revalidate();
      });
  },
  execute: function (options) {
    var file = options.input.file;
    var _this = this;
    phantasus.OpenDatasetTool
      .fileExtensionPrompt(file,
        function (readOptions) {
          if (!readOptions) {
            readOptions = {};
          }
          readOptions.interactive = true;
          var deferred = phantasus.DatasetUtil.read(file,
            readOptions);
          _this._read(options, deferred);
        });

  }, // prompt for metadata field name in dataset and in file
  _matchAppend: function (newDatasetMetadataNames,
                          currentDatasetMetadataNames, heatMap, callback) {
    var tool = {};
    tool.execute = function (options) {
      return options.input;
    };
    tool.toString = function () {
      return 'Select Fields';
    };
    tool.gui = function () {
      var items = [{
        name: 'current_dataset_annotation_name',
        options: currentDatasetMetadataNames,
        type: 'select',
        value: 'id',
        required: true
      }];
      items.push({
        name: 'new_dataset_annotation_name',
        type: 'select',
        value: 'id',
        options: newDatasetMetadataNames,
        required: true
      });
      return items;
    };
    phantasus.HeatMap.showTool(tool, heatMap, callback);
  },
  _matchOverlay: function (newDatasetColumnMetadataNames,
                           currentDatasetColumnMetadataNames, newDatasetRowMetadataNames,
                           currentDatasetRowMetadataNames, heatMap, callback) {
    var tool = {};
    tool.execute = function (options) {
      return options.input;
    };
    tool.toString = function () {
      return 'Select Fields';
    };
    tool.gui = function () {
      var items = [];
      items.push({
        name: 'current_dataset_column_annotation_name',
        options: currentDatasetColumnMetadataNames,
        type: 'select',
        value: 'id',
        required: true
      });
      items.push({
        name: 'new_dataset_column_annotation_name',
        type: 'select',
        value: 'id',
        options: newDatasetColumnMetadataNames,
        required: true
      });
      items.push({
        name: 'current_dataset_row_annotation_name',
        options: currentDatasetRowMetadataNames,
        type: 'select',
        value: 'id',
        required: true
      });
      items.push({
        name: 'new_dataset_row_annotation_name',
        type: 'select',
        value: 'id',
        options: newDatasetRowMetadataNames,
        required: true
      });
      return items;
    };
    phantasus.HeatMap.showTool(tool, heatMap, callback);
  }
};
