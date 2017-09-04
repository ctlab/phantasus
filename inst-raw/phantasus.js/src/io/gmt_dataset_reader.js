phantasus.GmtDatasetReader = function () {
};
phantasus.GmtDatasetReader.prototype = {
  getFormatName: function () {
    return 'gmt';
  },
  read: function (fileOrUrl, callback) {
    var name = phantasus.Util.getBaseFileName(phantasus.Util
      .getFileName(fileOrUrl));
    phantasus.ArrayBufferReader.getArrayBuffer(fileOrUrl, function (err, arrayBuffer) {
      if (err) {
        callback(err);
      } else {
        try {
          callback(null, phantasus.DatasetUtil.geneSetsToDataset(name,
            new phantasus.GmtReader()
              .read(new phantasus.ArrayBufferReader(
                new Uint8Array(arrayBuffer)))));
        }
        catch (x) {
          callback(x);
        }
      }
    });

  }
};
