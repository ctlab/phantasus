phantasus.GeoReader = function () {
};

phantasus.GeoReader.prototype = {
  read: function (name, callback) {
    var req = ocpu.call('loadGEO', { name: name }, function (session) {
      session.getObject(function (success) {
        //// console.log('phantasus.GeoReader.prototype.read ::', success);

        var filePath = phantasus.Util.getFilePath(session, JSON.parse(success)[0]);
        //// console.log(filePath);

        var r = new FileReader();

        r.onload = function (e) {
          var contents = e.target.result;
          //// console.log(contents);
          var ProtoBuf = dcodeIO.ProtoBuf;
          ProtoBuf.protoFromFile("./message.proto", function (error, success) {
            if (error) {
              alert(error);
              // console.log("GeoReader ::", "ProtoBuilder failed", error);
              return;
            }
            var builder = success,
              rexp = builder.build("rexp"),
              REXP = rexp.REXP,
              rclass = REXP.RClass;


            var res = REXP.decode(contents);

            var jsondata = phantasus.Util.getRexpData(res, rclass);
            //// console.log(jsondata);

            var datasets = [];
            for (var k = 0; k < Object.keys(jsondata).length; k++) {
              var dataset = phantasus.GeoReader.getDataset(session, Object.keys(jsondata)[k], jsondata[Object.keys(jsondata)[k]]);
              dataset.setESVariable('es_' + (k + 1).toString());
              datasets.push(dataset);
            }
            // console.log("resulting datasets", datasets);
            callback(null, datasets);
          });
        };

        phantasus.BlobFromPath.getFileObject(filePath, function (f) {
          r.readAsArrayBuffer(f);
        });
      })
    });
    req.fail(function () {
      callback(req.responseText);
      //// console.log('phantasus.GeoReader.prototype.read ::', req.responseText);
    });

  },
  _parse: function (text) {

  }
};


phantasus.GeoReader.getDataset = function (session, seriesName, jsondata) {
  var flatData = jsondata.data.values;
  var nrowData = jsondata.data.dim[0];
  var ncolData = jsondata.data.dim[1];
  var flatPdata = jsondata.pdata.values;
  //var participants = jsondata.participants.values;
  var annotation = jsondata.fdata.values;
  //// console.log(annotation);
  var id = jsondata.rownames.values;
  var metaNames = jsondata.colMetaNames.values;
  var rowMetaNames = jsondata.rowMetaNames.values;

  var matrix = [];
  for (var i = 0; i < nrowData; i++) {
    var curArray = new Float32Array(ncolData);
    for (var j = 0; j < ncolData; j++) {
      curArray[j] = flatData[i + j * nrowData];
    }
    matrix.push(curArray);
  }
  var dataset = new phantasus.Dataset({
    name: seriesName,
    rows: nrowData,
    columns: ncolData,
    array: matrix,
    dataType: 'Float32',
    esSession: session,
    isGEO: true
  });

  for (i = 0; i < metaNames.length; i++) {
    var curVec = dataset.getColumnMetadata().add(metaNames[i]);
    for (j = 0; j < ncolData; j++) {
      curVec.setValue(j, flatPdata[j + i * ncolData]);
    }
  }

  var rowIds = dataset.getRowMetadata().add('id');

  for (i = 0; i < rowMetaNames.length; i++) {
    curVec = dataset.getRowMetadata().add(rowMetaNames[i]);
    for (j = 0; j < nrowData; j++) {
      curVec.setValue(j, annotation[j + i * nrowData]);
      rowIds.setValue(j, id[j])
    }
  }
  phantasus.MetadataUtil.maybeConvertStrings(dataset.getRowMetadata(), 1);
  phantasus.MetadataUtil.maybeConvertStrings(dataset.getColumnMetadata(),
    1);

  //// console.log("returned dataset", dataset);
  return dataset;
};
