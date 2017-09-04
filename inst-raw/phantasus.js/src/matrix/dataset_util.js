/**
 * Static utilities for phantasus.DatasetInterface instances
 *
 * @class phantasus.DatasetUtil
 */
phantasus.DatasetUtil = function () {
};
phantasus.DatasetUtil.min = function (dataset, seriesIndex) {
  seriesIndex = seriesIndex || 0;
  var min = Number.MAX_VALUE;
  for (var i = 0, rows = dataset.getRowCount(); i < rows; i++) {
    for (var j = 0, columns = dataset.getColumnCount(); j < columns; j++) {
      var d = dataset.getValue(i, j, seriesIndex);
      if (isNaN(d)) {
        continue;
      }
      min = Math.min(min, d);
    }
  }
  return min;
};
phantasus.DatasetUtil.slicedView = function (dataset, rows, columns) {
  return new phantasus.SlicedDatasetView(dataset, rows, columns);
};
phantasus.DatasetUtil.transposedView = function (dataset) {
  return dataset instanceof phantasus.TransposedDatasetView ? dataset
    .getDataset() : new phantasus.TransposedDatasetView(dataset);
};
phantasus.DatasetUtil.max = function (dataset, seriesIndex) {
  seriesIndex = seriesIndex || 0;
  var max = -Number.MAX_VALUE;
  for (var i = 0, rows = dataset.getRowCount(); i < rows; i++) {
    for (var j = 0, columns = dataset.getColumnCount(); j < columns; j++) {
      var d = dataset.getValue(i, j, seriesIndex);
      if (isNaN(d)) {
        continue;
      }
      max = Math.max(max, d);
    }
  }
  return max;
};

phantasus.DatasetUtil.getDatasetReader = function (ext, options) {
  if (options == null) {
    options = {};
  }
  var datasetReader;
  if (ext === 'maf') {
    datasetReader = new phantasus.MafFileReader();
    if (options && options.mafGeneFilter) {
      datasetReader.setGeneFilter(options.mafGeneFilter);
    }
  } else if (ext === 'gct') {
    datasetReader = new phantasus.GctReader();
    // datasetReader = new phantasus.StreamingGctReader();
  } else if (ext === 'gmt') {
    datasetReader = new phantasus.GmtDatasetReader();
  } else if (ext === 'xlsx' || ext === 'xls') {
    datasetReader = options.interactive ? new phantasus.Array2dReaderInteractive() : new phantasus.XlsxDatasetReader();
  } else if (ext === 'segtab' || ext === 'seg') {
    datasetReader = new phantasus.SegTabReader();
    if (options && options.regions) {
      datasetReader.setRegions(options.regions);
    }
  } else if (ext === 'txt' || ext === 'tsv' || ext === 'csv') {
    datasetReader = options.interactive ? new phantasus.Array2dReaderInteractive() : new phantasus.TxtReader();
  } else if (ext === 'json') {
    datasetReader = new phantasus.JsonDatasetReader();
  } else {
    datasetReader = new phantasus.GctReader();
  }
  return datasetReader;
}
;

phantasus.DatasetUtil.readDatasetArray = function (datasets) {
  var retDef = $.Deferred();
  var loadedDatasets = [];
  var promises = [];
  _.each(datasets, function (url, i) {
    var p = phantasus.DatasetUtil.read(url);
    p.index = i;
    p.done(function (dataset) {
      loadedDatasets[this.index] = dataset;
    });
    p.fail(function (err) {
      var message = ['Error opening ' + phantasus.Util
        .getFileName(url) + '.'];
      if (err.message) {
        message.push('<br />Cause: ');
        message.push(err.message);
      }
      retDef.reject(message.join(''));

    });
    promises.push(p);
  });
  if (promises.length === 0) {
    retDef.reject('No datasets specified.');
  }

  $.when
    .apply($, promises)
    .then(
      function () {
        retDef.resolve(phantasus.DatasetUtil.join(loadedDatasets, 'id'));
      });
  return retDef;
};
/**
 * Annotate a dataset from external file or text.
 *
 * @param options.annotations -
 *            Array of file, datasetField, and fileField.
 * @param options.isColumns -
 *            Whether to annotate columns
 * @return A jQuery Deferred object that resolves to an array of functions to
 *         execute with a dataset parameter.
 */
phantasus.DatasetUtil.annotate = function (options) {
  var retDef = $.Deferred();
  var promises = [];
  var functions = [];
  var isColumns = options.isColumns;
  _.each(options.annotations, function (ann, annotationIndex) {
    if (phantasus.Util.isArray(ann.file)) { // already parsed text
      functions[annotationIndex] = function (dataset) {
        new phantasus.OpenFileTool().annotate(ann.file, dataset,
          isColumns, null, ann.datasetField, ann.fileField,
          ann.include);
      };
    } else {
      var result = phantasus.Util.readLines(ann.file);
      var fileName = phantasus.Util.getFileName(ann.file);
      var deferred = $.Deferred();
      promises.push(deferred);
      result.fail(function (message) {
        deferred.reject(message);
      });
      result.done(function (lines) {
        if (phantasus.Util.endsWith(fileName, '.gmt')) {
          var sets = new phantasus.GmtReader().parseLines(lines);
          functions[annotationIndex] = function (dataset) {
            new phantasus.OpenFileTool().annotate(null, dataset,
              isColumns, sets, ann.datasetField,
              ann.fileField);
          };
          deferred.resolve();
        } else if (phantasus.Util.endsWith(fileName, '.cls')) {
          functions[annotationIndex] = function (dataset) {
            new phantasus.OpenFileTool().annotateCls(null, dataset,
              fileName, isColumns, lines);
          };
          deferred.resolve();
        } else {
          functions[annotationIndex] = function (dataset) {
            new phantasus.OpenFileTool().annotate(lines, dataset,
              isColumns, null, ann.datasetField,
              ann.fileField, ann.include);
          };
          deferred.resolve();
        }
      });
    }
  });
  $.when.apply($, promises).then(function () {
    retDef.resolve(functions);
  });
  return retDef;
};
/**
 * Reads a dataset at the specified URL or file
 * @param file
 *            a File or URL
 * @return A promise that resolves to phantasus.DatasetInterface
 */
phantasus.DatasetUtil.read = function (fileOrUrl, options) {
  if (fileOrUrl == null) {
    throw 'File is null';
  }
  if (options == null) {
    options = {};
  }
  var isFile = fileOrUrl instanceof File;
  var isString = phantasus.Util.isString(fileOrUrl);
  var ext = options.extension ? options.extension : phantasus.Util.getExtension(phantasus.Util.getFileName(fileOrUrl));
  var datasetReader;
  var str = fileOrUrl.toString();

  var isGSE = isString && (fileOrUrl.substring(0, 3) === 'GSE' || fileOrUrl.substring(0, 3) === 'GDS');

  if (isGSE) {
    datasetReader = new phantasus.GeoReader({type: fileOrUrl.substring(0, 3)});
  }
  else if (ext === '' && str != null && str.indexOf('blob:') === 0) {
    datasetReader = new phantasus.TxtReader(); // copy from clipboard
  } else {
    datasetReader = phantasus.DatasetUtil.getDatasetReader(ext, options);
  }
  if (isString || isFile) { // URL or file
    var deferred = $.Deferred();
    // override toString so can determine file name
    if (options.background) {
      var path = phantasus.Util.getScriptPath();
      var blob = new Blob(
        ['self.onmessage = function(e) {'
        + 'importScripts(e.data.path);'
        + 'var ext = phantasus.Util.getExtension(phantasus.Util'
        + '.getFileName(e.data.fileOrUrl));'
        + 'var datasetReader = phantasus.DatasetUtil.getDatasetReader(ext,'
        + '	e.data.options);'
        + 'datasetReader.read(e.data.fileOrUrl, function(err,dataset) {'
        + '	self.postMessage(dataset);' + '	});' + '}']);

      var blobURL = window.URL.createObjectURL(blob);
      var worker = new Worker(blobURL);
      worker.addEventListener('message', function (e) {
        deferred.resolve(phantasus.Dataset.fromJSON(e.data));
        window.URL.revokeObjectURL(blobURL);
      }, false);
      // start the worker
      worker.postMessage({
        path: path,
        fileOrUrl: fileOrUrl,
        options: options
      });

    } else {
      datasetReader.read(fileOrUrl, function (err, dataset) {
        if (err) {
          deferred.reject(err);
        } else {

          // console.log(dataset);
          // console.log('ready to resolve with', dataset);
          deferred.resolve(dataset);
          if (!isGSE) {
            phantasus.DatasetUtil.toESSessionPromise({dataset: dataset, isGEO: isGSE});
          }
        }
      });

    }
    var pr = deferred.promise();
    pr.toString = function () {
      return '' + fileOrUrl;
    };
    return pr;
  } else if (typeof fileOrUrl.done === 'function') { // assume it's a
    // deferred
    return fileOrUrl;
  } else { // it's already a dataset?
    var deferred = $.Deferred();
    if (fileOrUrl.getRowCount) {
      deferred.resolve(fileOrUrl);
    } else {
      deferred.resolve(phantasus.Dataset.fromJSON(fileOrUrl));
    }
    return deferred.promise();
  }

};

/**
 * @param dataset
 *            The dataset to convert to an array
 * @param options.columns
 *            An array of column indices to include from the dataset
 * @param options.columnFields
 *            An array of field names to use in the returned objects that
 *            correspond to the column indices in the dataset
 * @param options.metadataFields
 *            An array of row metadata fields to include from the dataset
 *
 */
phantasus.DatasetUtil.toObjectArray = function (dataset, options) {
  var columns = options.columns || [0];
  var columnFields = options.columnFields || ['value'];
  if (columnFields.length !== columns.length) {
    throw 'columns.length !== columnFields.length';
  }
  var metadataFields = options.metadataFields;
  // grab all of the headers and filter the meta data vectors in the dataset
  // down
  // to the ones specified in metaFields. If metaFields is not passed, take
  // all metadata
  var rowMetadata = dataset.getRowMetadata();
  if (!metadataFields) {
    metadataFields = phantasus.MetadataUtil.getMetadataNames(rowMetadata);
  }
  var vectors = phantasus.MetadataUtil.getVectors(rowMetadata, metadataFields);
  // build an object that contains the matrix values for the given columns
  // along
  // with any metadata
  var array = [];
  for (var i = 0; i < dataset.getRowCount(); i++) {
    var obj = {};
    for (var j = 0; j < columns.length; j++) {
      obj[columnFields[j]] = dataset.getValue(i, columns[j]);
    }
    for (var j = 0; j < vectors.length; j++) {
      obj[vectors[j].getName()] = vectors[j].getValue(i);
    }
    array.push(obj);
  }
  return array;
};
phantasus.DatasetUtil.fixL1K = function (dataset) {
  var names = {
    'cell_id': 'Cell Line',
    'pert_idose': 'Dose (\u00B5M)',
    'pert_iname': 'Name',
    'pert_itime': 'Time (hr)',
    'distil_ss': 'Signature Strength',
    'pert_type': 'Type',
    'cell_lineage': 'Lineage',
    'cell_histology': 'Histology',
    'cell_type': 'Cell Type'
  };
  var fixNames = function (metadata) {
    for (var i = 0, count = metadata.getMetadataCount(); i < count; i++) {
      var v = metadata.get(i);
      var name = v.getName();
      var mapped = names[name];
      if (mapped) {
        v.setName(mapped);
      }
    }
  };
  fixNames(dataset.getRowMetadata());
  fixNames(dataset.getColumnMetadata());
  var fix666 = function (metadata) {
    for (var i = 0, count = metadata.getMetadataCount(); i < count; i++) {
      var v = metadata.get(i);
      if (v.getName() == 'Dose (\u00B5M)') { // convert to number
        for (var j = 0, size = v.size(); j < size; j++) {
          var value = v.getValue(j);
          if (value != null) {
            v.setValue(j, parseFloat(value));
          }
        }
      }
      var isNumber = false;
      for (var j = 0, size = v.size(); j < size; j++) {
        var value = v.getValue(j);
        if (value != null) {
          isNumber = _.isNumber(value);
          break;
        }
      }
      var newValue = isNumber || v.getName() == 'Dose (\u00B5M)' ? 0 : '';
      for (var j = 0, size = v.size(); j < size; j++) {
        var value = v.getValue(j);
        if (value != null && value == '-666') {
          v.setValue(j, newValue);
        }
      }
    }
  };
  fix666(dataset.getRowMetadata());
  fix666(dataset.getColumnMetadata());
  var fixCommas = function (metadata) {
    var regex = /(,)([^ ])/g;
    _.each(['Lineage', 'Histology'], function (name) {
      var v = metadata.getByName(name);
      if (v != null) {
        for (var i = 0, size = v.size(); i < size; i++) {
          var val = v.getValue(i);
          if (val) {
            v.setValue(i, val.replace(regex, ', $2'));
          }
        }
      }
    });
  };
  fixCommas(dataset.getRowMetadata());
  fixCommas(dataset.getColumnMetadata());
};
phantasus.DatasetUtil.geneSetsToDataset = function (name, sets) {
  var uniqueIds = new phantasus.Map();
  for (var i = 0, length = sets.length; i < length; i++) {
    var ids = sets[i].ids;
    for (var j = 0, nIds = ids.length; j < nIds; j++) {
      uniqueIds.set(ids[j], 1);
    }
  }
  var uniqueIdsArray = uniqueIds.keys();
  var dataset = new phantasus.Dataset({
    name: name,
    rows: uniqueIdsArray.length,
    columns: sets.length
  });
  var columnIds = dataset.getColumnMetadata().add('id');
  for (var i = 0, length = sets.length; i < length; i++) {
    columnIds.setValue(i, sets[i].name);
  }
  var rowIds = dataset.getRowMetadata().add('id');
  for (var i = 0, size = uniqueIdsArray.length; i < size; i++) {
    rowIds.setValue(i, uniqueIdsArray[i]);
  }
  var rowIdToIndex = phantasus.VectorUtil.createValueToIndexMap(rowIds);
  for (var i = 0, length = sets.length; i < length; i++) {
    var ids = sets[i].ids;
    for (var j = 0, nIds = ids.length; j < nIds; j++) {
      dataset.setValue(rowIdToIndex.get(ids[j]), i, 1);
    }
  }
  return dataset;
};
phantasus.DatasetUtil.DATASET_FILE_FORMATS = '<a target="_blank" href="https://clue.io/help#datasets">GCT 1.3</a>, '
  + '<a target="_blank" href="http://www.broadinstitute.org/cancer/software/genepattern/gp_guides/file-formats/sections/gct">GCT 1.2</a>, '
  + '<a target="_blank" href="https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification">MAF</a>, '
  + '<a target="_blank" href="http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29">GMT</a>, '
  + ' a tab-delimited text file, or an Excel spreadsheet';
phantasus.DatasetUtil.SESSION_FILE_FORMAT = 'a saved phantasus session';

phantasus.DatasetUtil.DATASET_AND_SESSION_FILE_FORMATS = '<a target="_blank"' +
  ' href="https://clue.io/help#datasets">GCT 1.3</a>, '
  + '<a target="_blank" href="http://www.broadinstitute.org/cancer/software/genepattern/gp_guides/file-formats/sections/gct">GCT 1.2</a>, '
  + '<a target="_blank" href="https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification">MAF</a>, '
  + '<a target="_blank" href="http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29">GMT</a>, '
  + ' a tab-delimited text file, an Excel spreadsheet, or a saved phantasus session';
phantasus.DatasetUtil.BASIC_DATASET_FILE_FORMATS = '<a target="_blank" href="https://clue.io/help#datasets">GCT 1.3</a>, '
  + '<a target="_blank" href="http://www.broadinstitute.org/cancer/software/genepattern/gp_guides/file-formats/sections/gct">GCT 1.2</a>, '
  + ' or a tab-delimited text file';
phantasus.DatasetUtil.GCT_FILE_FORMAT = '<a target="_blank" href="https://clue.io/help#datasets">GCT 1.3</a>';
phantasus.DatasetUtil.ANNOTATION_FILE_FORMATS = 'an xlsx file, tab-delimited text file, or a <a target="_blank" href="http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29">GMT file</a>';
phantasus.DatasetUtil.DENDROGRAM_FILE_FORMATS = 'a <a href="http://en.wikipedia.org/wiki/Newick_format" target="_blank">Newick</a> file';
phantasus.DatasetUtil.OPEN_FILE_FORMATS = '<a target="_blank" href="https://clue.io/help#datasets">GCT 1.3</a>, '
  + '<a target="_blank" href="http://www.broadinstitute.org/cancer/software/genepattern/gp_guides/file-formats/sections/gct">GCT 1.2</a>, '
  + '<a target="_blank" href="https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification">MAF</a>, '
  + '<a target="_blank" href="http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29">GMT</a>, '
  + ' a tab-delimited text file, or a <a href="http://en.wikipedia.org/wiki/Newick_format" target="_blank">Newick</a> file';
phantasus.DatasetUtil.getRootDataset = function (dataset) {
  while (dataset.getDataset) {
    dataset = dataset.getDataset();
  }
  return dataset;
};

phantasus.DatasetUtil.getSeriesIndex = function (dataset, name) {
  for (var i = 0, nseries = dataset.getSeriesCount(); i < nseries; i++) {
    if (name === dataset.getName(i)) {
      return i;
    }
  }
  return -1;
};
phantasus.DatasetUtil.getSeriesNames = function (dataset) {
  var names = [];
  for (var i = 0, nseries = dataset.getSeriesCount(); i < nseries; i++) {
    names.push(dataset.getName(i));
  }
  // names.sort(function (a, b) {
  // 	a = a.toLowerCase();
  // 	b = b.toLowerCase();
  // 	return (a < b ? -1 : (a === b ? 0 : 1));
  // });
  return names;
};

/**
 * Search dataset values.
 *
 * @param options.dataset
 *      The dataset
 * @param options.text
 *            Search text
 * @param options.defaultMatchMode
 *            'exact' or 'contains'
 * @param options.matchAllPredicates Whether to match all predicates
 * @return Set of matching indices.
 *
 */
phantasus.DatasetUtil.searchValues = function (options) {
  if (text === '') {
    return;
  }
  var dataset = options.dataset;
  var text = options.text;
  var tokens = phantasus.Util.getAutocompleteTokens(text);
  if (tokens.length == 0) {
    return;
  }
  var predicates = phantasus.Util.createSearchPredicates({
    tokens: tokens,
    defaultMatchMode: options.defaultMatchMode
  });
  var matchAllPredicates = options.matchAllPredicates === true;
  var npredicates = predicates.length;
  var viewIndices = new phantasus.Set();

  function isMatch(object, toObject, predicate) {
    if (object != null) {
      if (toObject) {
        var filterColumnName = predicate.getField();
        if (filterColumnName != null) {
          var value = object[filterColumnName];
          return predicate.accept(value);
        } else { // try all fields
          for (var name in object) {
            var value = object[name];
            return predicate.accept(value);
          }
        }
      } else {
        var filterColumnName = predicate.getField();
        if (filterColumnName == null || filterColumnName === dataset.getName(k)) {
          return predicate.accept(object);

        }
      }
    }
  }

  for (var i = 0, nrows = dataset.getRowCount(); i < nrows; i++) {
    for (var j = 0, ncols = dataset.getColumnCount(); j < ncols; j++) {
      var matches = false;
      itemSearch:
        if (matchAllPredicates) {
          matches = true;
          for (var p = 0; p < npredicates; p++) {
            var predicate = predicates[p];
            var pmatch = false;
            for (var k = 0, nseries = dataset.getSeriesCount(); k < nseries; k++) {
              var element = dataset.getValue(i, j, k);
              var isObject = element != null && element.toObject != null;
              if (isObject) {
                element = element.toObject();
              }
              if (isMatch(element, isObject, predicate)) {
                pmatch = true;
                break;
              }
            }
            if (!pmatch) {
              matches = false;
              break itemSearch;
            }
          }
        } else {
          for (var k = 0, nseries = dataset.getSeriesCount(); k < nseries; k++) {
            var element = dataset.getValue(i, j, k);
            var isObject = element != null && element.toObject != null;
            if (isObject) {
              element = element.toObject();
            }
            for (var p = 0; p < npredicates; p++) {
              var predicate = predicates[p];
              if (isMatch(element, isObject, predicate)) {
                matches = true;
                break itemSearch;
              }
            }
          }
        }

      if (matches) {
        viewIndices
          .add(new phantasus.Identifier(
            [i, j]));
      }
    }
  }
  return viewIndices;

};

/**
 * Search dataset values.
 */
phantasus.DatasetUtil.autocompleteValues = function (dataset) {
  return function (tokens, cb) {

    var token = tokens != null && tokens.length > 0 ? tokens[tokens.selectionStartIndex]
      : '';
    token = $.trim(token);
    var seriesIndices = [];
    for (var i = 0, nrows = dataset.getRowCount(); i < nrows; i++) {
      for (var k = 0, nseries = dataset.getSeriesCount(); k < nseries; k++) {
        if (dataset.getDataType(i, k) === 'Number') {
          seriesIndices.push([i, k]);
        }
      }
    }
    if (seriesIndices.length === 0) {
      return cb();
    }
    var _val; // first non-null value
    elementSearch: for (var k = 0, nseries = seriesIndices.length; k < nseries; k++) {
      var pair = seriesIndices[k];
      for (var j = 0, ncols = dataset.getColumnCount(); j < ncols; j++) {
        var element = dataset.getValue(pair[0], j, pair[1]);
        if (element != null && element.toObject) {
          _val = element.toObject();
          break elementSearch;
        }
      }
    }
    var matches = [];
    var fields = _val == null ? [] : _.keys(_val);
    if (token === '') {
      fields.sort(function (a, b) {
        return (a === b ? 0 : (a < b ? -1 : 1));
      });
      fields.forEach(function (field) {
        matches.push({
          value: field + ':',
          label: '<span style="font-weight:300;">' + field
          + ':</span>',
          show: true
        });
      });
      return cb(matches);
    }

    var field = null;
    var semi = token.indexOf(':');
    if (semi > 0) { // field search?
      if (token.charCodeAt(semi - 1) !== 92) { // \:
        var possibleField = $.trim(token.substring(0, semi));
        if (possibleField.length > 0 && possibleField[0] === '"'
          && possibleField[token.length - 1] === '"') {
          possibleField = possibleField.substring(1,
            possibleField.length - 1);
        }
        var index = fields.indexOf(possibleField);
        if (index !== -1) {
          token = $.trim(token.substring(semi + 1));
          field = possibleField;
        }
      }

    }

    var set = new phantasus.Set();
    // regex used to determine if a string starts with substring `q`
    var regex = new RegExp('^' + phantasus.Util.escapeRegex(token), 'i');
    // iterate through the pool of strings and for any string that
    // contains the substring `q`, add it to the `matches` array
    var max = 10;

    loop: for (var k = 0, nseries = seriesIndices.length; k < nseries; k++) {
      var pair = seriesIndices[k];
      for (var j = 0, ncols = dataset.getColumnCount(); j < ncols; j++) {
        var element = dataset.getValue(pair[0], j, pair[1]);
        if (element && element.toObject) {
          var object = element.toObject();
          if (field !== null) {
            var val = object[field];
            if (val != null) {
              var id = new phantasus.Identifier([val, field]);
              if (!set.has(id) && regex.test(val)) {
                set.add(id);
                if (set.size() === max) {
                  break loop;
                }
              }
            }
          } else { // search all fields
            for (var name in object) {
              var val = object[name];
              var id = new phantasus.Identifier([val, name]);
              if (!set.has(id) && regex.test(val)) {
                set.add(id);
                if (set.size() === max) {
                  break loop;
                }
              }
            }
          }

        }
      }
    }
    set.forEach(function (id) {
      var array = id.getArray();
      var field = array[1];
      var val = array[0];
      matches.push({
        value: field + ':' + val,
        label: '<span style="font-weight:300;">' + field + ':</span>'
        + '<span style="font-weight:900;">' + val + '</span>'
      });

    });
    if (field == null) {
      fields.forEach(function (field) {
        if (regex.test(field)) {
          matches.push({
            value: field + ':',
            label: '<span style="font-weight:300;">' + field
            + ':</span>',
            show: true
          });
        }
      });
    }
    cb(matches);
  };

};
// phantasus.DatasetUtil.toJSON = function(dataset) {
// var json = [];
// json.push('{');
// json.push('"name":"' + dataset.getName() + '", ');
// json.push('"v":['); // row major 2d array
// for (var i = 0, nrows = dataset.getRowCount(); i < nrows; i++) {
// if (i > 0) {
// json.push(',\n');
// }
// json.push('[');
// for (var j = 0, ncols = dataset.getColumnCount(); j < ncols; j++) {
// if (j > 0) {
// json.push(',');
// }
// json.push(JSON.stringify(dataset.getValue(i, j)));
// }
// json.push(']');
// }
// json.push(']'); // end v
// var metadatatoJSON = function(model) {
// json.push('[');
// for (var i = 0, count = model.getMetadataCount(); i < count; i++) {
// var v = model.get(i);
// if (i > 0) {
// json.push(',\n');
// }
// json.push('{');
// json.push('"id":"' + v.getName() + '"');
// json.push(', "v":[');
// for (var j = 0, nitems = v.size(); j < nitems; j++) {
// if (j > 0) {
// json.push(',');
// }
// json.push(JSON.stringify(v.getValue(j)));
// }
// json.push(']'); // end v array
// json.push('}');
// }
// json.push(']');
// };
// json.push(', "cols":');
// metadatatoJSON(dataset.getColumnMetadata());
// json.push(', "rows":');
// metadatatoJSON(dataset.getRowMetadata());
// json.push('}'); // end json object
// return json.join('');
// };
phantasus.DatasetUtil.fill = function (dataset, value, seriesIndex) {
  seriesIndex = seriesIndex || 0;
  for (var i = 0, nrows = dataset.getRowCount(), ncols = dataset
    .getColumnCount(); i < nrows; i++) {
    for (var j = 0; j < ncols; j++) {
      dataset.setValue(i, j, value, seriesIndex);
    }
  }
};

/**
 * Add an additional series to a dataset from another dataset.
 * @param options.dataset The dataset to add a series to
 * @param options.newDataset The dataset that is used as the source for the overlay
 * @param options.rowAnnotationName dataset row annotation name to use for matching
 * @param options.columnAnnotationName dataset column annotation name to use for matching
 * @param options.newRowAnnotationName newDataset row annotation name to use for matching
 * @param options.newColumnAnnotationName newDataset column annotation name to use for matching
 *
 */
phantasus.DatasetUtil.overlay = function (options) {
  var dataset = options.dataset;
  var newDataset = options.newDataset;
  var current_dataset_row_annotation_name = options.rowAnnotationName;
  var current_dataset_column_annotation_name = options.columnAnnotationName;
  var new_dataset_row_annotation_name = options.newRowAnnotationName;
  var new_dataset_column_annotation_name = options.newColumnAnnotationName;

  var rowValueToIndexMap = phantasus.VectorUtil
    .createValueToIndexMap(dataset
      .getRowMetadata()
      .getByName(
        current_dataset_row_annotation_name));
  var columnValueToIndexMap = phantasus.VectorUtil
    .createValueToIndexMap(dataset
      .getColumnMetadata()
      .getByName(
        current_dataset_column_annotation_name));
  var seriesIndex = dataset
    .addSeries({
      name: newDataset
        .getName(),
      dataType: newDataset.getDataType(0)
    });

  var rowVector = newDataset
    .getRowMetadata()
    .getByName(
      new_dataset_row_annotation_name);
  var rowIndices = [];
  var newDatasetRowIndicesSubset = [];
  for (var i = 0, size = rowVector
    .size(); i < size; i++) {
    var index = rowValueToIndexMap
      .get(rowVector
        .getValue(i));
    if (index !== undefined) {
      rowIndices.push(index);
      newDatasetRowIndicesSubset
        .push(i);
    }
  }

  var columnVector = newDataset
    .getColumnMetadata()
    .getByName(
      new_dataset_column_annotation_name);
  var columnIndices = [];
  var newDatasetColumnIndicesSubset = [];
  for (var i = 0, size = columnVector
    .size(); i < size; i++) {
    var index = columnValueToIndexMap
      .get(columnVector
        .getValue(i));
    if (index !== undefined) {
      columnIndices.push(index);
      newDatasetColumnIndicesSubset
        .push(i);
    }
  }
  newDataset = new phantasus.SlicedDatasetView(
    newDataset,
    newDatasetRowIndicesSubset,
    newDatasetColumnIndicesSubset);
  for (var i = 0, nrows = newDataset
    .getRowCount(); i < nrows; i++) {
    for (var j = 0, ncols = newDataset
      .getColumnCount(); j < ncols; j++) {
      dataset.setValue(
        rowIndices[i],
        columnIndices[j],
        newDataset
          .getValue(
            i,
            j),
        seriesIndex);

    }
  }
};
/**
 * Joins datasets by appending rows.
 * @param datasets
 * @param field
 * @return {phantasus.AbstractDataset} The joined dataset.
 */
phantasus.DatasetUtil.join = function (datasets, field) {
  if (datasets.length === 0) {
    throw 'No datasets';
  }
  if (datasets.length === 1) {
    var name = datasets[0].getName();
    var sourceVector = datasets[0].getRowMetadata().add('Source');
    for (var i = 0, size = sourceVector.size(); i < size; i++) {
      sourceVector.setValue(i, name);
    }
    return datasets[0];
  }
  // take union of all ids
  var ids = new phantasus.Set();
  for (var i = 0; i < datasets.length; i++) {
    var idVector = datasets[i].getColumnMetadata().getByName(field);
    for (var j = 0, size = idVector.size(); j < size; j++) {
      ids.add(idVector.getValue(j));
    }
  }
  var dummyDataset = new phantasus.Dataset({
    rows: 0,
    columns: ids.size(),
    name: datasets[0].getName()
  });
  var dummyIdVector = dummyDataset.getColumnMetadata().add(field);
  var counter = 0;
  ids.forEach(function (id) {
    dummyIdVector.setValue(counter++, id);
  });

  var dataset = new phantasus.JoinedDataset(
    dummyDataset, datasets[0], field,
    field);
  for (var i = 1; i < datasets.length; i++) {
    dataset = new phantasus.JoinedDataset(dataset,
      datasets[i], field, field);
  }
  return dataset;
};
phantasus.DatasetUtil.shallowCopy = function (dataset) {
  // make a shallow copy of the dataset, metadata is immutable via the UI
  var rowMetadataModel = phantasus.MetadataUtil.shallowCopy(dataset
    .getRowMetadata());
  var columnMetadataModel = phantasus.MetadataUtil.shallowCopy(dataset
    .getColumnMetadata());
  dataset.getRowMetadata = function () {
    return rowMetadataModel;
  };
  dataset.getColumnMetadata = function () {
    return columnMetadataModel;
  };
  return dataset;
};

phantasus.DatasetUtil.copy = function (dataset) {
  var newDataset = new phantasus.Dataset({
    name: dataset.getName(),
    rows: dataset.getRowCount(),
    columns: dataset.getColumnCount(),
    dataType: dataset.getDataType(0)
  });
  for (var seriesIndex = 0,
         nseries = dataset.getSeriesCount(); seriesIndex < nseries; seriesIndex++) {
    if (seriesIndex > 0) {
      newDataset.addSeries({
        name: dataset.getName(seriesIndex),
        rows: dataset.getRowCount(),
        columns: dataset.getColumnCount(),
        dataType: dataset.getDataType(seriesIndex)
      });
    }
    for (var i = 0, nrows = dataset.getRowCount(), ncols = dataset
      .getColumnCount(); i < nrows; i++) {
      for (var j = 0; j < ncols; j++) {
        newDataset.setValue(i, j, dataset.getValue(i, j, seriesIndex),
          seriesIndex);
      }
    }
  }
  var rowMetadataModel = phantasus.MetadataUtil.shallowCopy(dataset
    .getRowMetadata());
  var columnMetadataModel = phantasus.MetadataUtil.shallowCopy(dataset
    .getColumnMetadata());
  newDataset.getRowMetadata = function () {
    return rowMetadataModel;
  };
  newDataset.getColumnMetadata = function () {
    return columnMetadataModel;
  };
  if (dataset.getESSession()) {
    newDataset.setESSession(dataset.getESSession());
  }
  return newDataset;
};
phantasus.DatasetUtil.toString = function (dataset, value, seriesIndex) {
  seriesIndex = seriesIndex || 0;
  var s = [];
  for (var i = 0, nrows = dataset.getRowCount(), ncols = dataset
    .getColumnCount(); i < nrows; i++) {
    for (var j = 0; j < ncols; j++) {
      if (j > 0) {
        s.push(', ');
      }
      s.push(phantasus.Util.nf(dataset.getValue(i, j, seriesIndex)));
    }
    s.push('\n');
  }
  return s.join('');
};
phantasus.DatasetUtil.getNonEmptyRows = function (dataset) {
  var rowsToKeep = [];
  for (var i = 0, nrows = dataset.getRowCount(); i < nrows; i++) {
    var keep = false;
    for (var j = 0, ncols = dataset.getColumnCount(); j < ncols; j++) {
      var value = dataset.getValue(i, j);
      if (!isNaN(value)) {
        keep = true;
        break;
      }
    }
    if (keep) {
      rowsToKeep.push(i);
    }
  }
  return rowsToKeep;
};
phantasus.DatasetUtil.getContentArray = function (dataset) {

  var array = [];
  var nr = dataset.rows;
  var nc = dataset.columns;
  //// console.log("getContentArray ::", "dataset:", dataset, "rows:", nr, "columns:", nc);

  for (var i = 0; i < nc; i++) {
    for (var j = 0; j < nr; j++) {
      array.push(dataset.getValue(j, i));
    }
  }
  //// console.log("getContentArray ::", array);
  return array;
};
phantasus.DatasetUtil.getMetadataArray = function (dataset) {
  var pDataArray = [];
  var labelDescription = [];
  //// console.log("phantasus.DatasetUtil.getMetadataArray ::", dataset);
  var columnMeta = dataset.getColumnMetadata();
  var features = columnMeta.getMetadataCount();
  var participants = dataset.getColumnCount();


  for (var j = 0; j < features; j++) {
    var vecJ = columnMeta.get(j);
    for (var l = 0; l < participants; l++) {
      pDataArray.push({
        strval: vecJ.getValue(l) ? vecJ.getValue(l).toString() : "",
        isNA: false
      });
    }
    labelDescription.push({
      strval: vecJ.getName(),
      isNA: false
    });

  }

  var rowMeta = dataset.getRowMetadata();
  var fDataArray = [];
  var varLabels = [];
  for (var j = 0; j < rowMeta.getMetadataCount(); j++) {
    var vecJ = rowMeta.get(j);
    for (var l = 0; l < dataset.getRowCount(); l++){
      fDataArray.push({
        strval: vecJ.getValue(l) ? vecJ.getValue(l).toString() : "",
        isNA: false
      });
    }
    varLabels.push({
      strval: vecJ.getName(),
      isNA: false
    });
  }

  return {
    pdata: pDataArray,
    varLabels: labelDescription,
    fdata: fDataArray,
    fvarLabels: varLabels
  };
};

phantasus.DatasetUtil.toESSessionPromise = function (options) {
  var dataset = options.dataset ? options.dataset : options;

  //// console.log("ENTERED TO_ESSESSION_PROMISE", dataset, options);
  //var copiedDataset = phantasus.DatasetUtil.copy(dataset);
  //// console.log("EsSessionPromise ::", "after copying", dataset);
  while (dataset.dataset) {
    dataset = dataset.dataset;
  }
  dataset.setESSession(new Promise(function (resolve, reject) {
    //// console.log("phantasus.DatasetUtil.toESSessionPromise ::", dataset, dataset instanceof phantasus.Dataset, dataset instanceof phantasus.SlicedDatasetView);
    /*		if (dataset.dataset) {
     //// console.log("phantasus.DatasetUtil.toESSessionPromise ::", "dataset in instanceof phantasus.SlicedDatasetView", "go deeper");
     phantasus.DatasetUtil.toESSessionPromise(dataset.dataset);
     }*/
    //// console.log("before going further", options);
    if (options.isGEO) {
      //// console.log("toESSession::", "resolving as geo dataset");
      resolve(dataset.getESSession());
      return;
    }

    var array = phantasus.DatasetUtil.getContentArray(dataset);
    var meta = phantasus.DatasetUtil.getMetadataArray(dataset);

    // console.log(array, array.length);
    var messageJSON = {
      rclass: "LIST",
      rexpValue: [{
        rclass: "REAL",
        realValue: array,
        attrName: "dim",
        attrValue: {
          rclass: "INTEGER",
          intValue: [dataset.getRowCount(), dataset.getColumnCount()]
        }
      }, {
        rclass: "STRING",
        stringValue: meta.pdata,
        attrName: "dim",
        attrValue: {
          rclass: "INTEGER",
          intValue: [dataset.getColumnCount(), meta.varLabels.length]
        }
      }, {
        rclass: "STRING",
        stringValue: meta.varLabels
      }, {
        rclass: "STRING",
        stringValue: meta.fdata,
        attrName: "dim",
        attrValue: {
          rclass: "INTEGER",
          intValue: [dataset.getRowCount(), meta.fvarLabels.length]
        }
      }, {
        rclass: "STRING",
        stringValue: meta.fvarLabels
      }],
      attrName: "names",
      attrValue: {
        rclass: "STRING",
        stringValue: [{
          strval: "data",
          isNA: false
        }, {
          strval: "pData",
          isNA: false
        }, {
          strval: "varLabels",
          isNA: false
        }, {
          strval: "fData",
          isNA: false
        }, {
          strval: "fvarLabels",
          isNA: false
        }]
      }
    };
    var ProtoBuf = dcodeIO.ProtoBuf;
    ProtoBuf.protoFromFile('./message.proto', function (error, success) {
      if (error) {
        alert(error);
        // console.log("ExpressionSetCreation :: ", "ProtoBuilder failed", error);
        return;
      }
      //// console.log("phantasus.DatasetUtil.toESSessionPromise ::", "protobuilder error", error);
      //// console.log("phantasus.DatasetUtil.toESSessionPromise ::", "protobuilder success", success);
      var builder = success,
        rexp = builder.build('rexp'),
        REXP = rexp.REXP;

      var proto = new REXP(messageJSON);
      var req = ocpu.call('createES', proto, function (session) {
        //// console.log("phantasus.DatasetUtil.toESSessionPromise ::", "from successful request", session);
        dataset.setESVariable('es');
        resolve(session);
      }, true);

      req.fail(function () {
        reject(req.responseText);
      });
    });
  }));
  /*var blob = new Blob([new Uint8Array((new REXP(messageJSON)).toArrayBuffer())], {type: "application/octet-stream"});
   saveAs(blob, "test1.bin");*/
};