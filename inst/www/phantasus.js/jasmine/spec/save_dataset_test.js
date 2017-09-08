describe('save_dataset_test', function () {

  it('save_dataset_to_gct', function () {
    var dataset;

    var promises = [];
    promises.push(phantasus.DatasetUtil.read('test_files/all_aml_train.gct').done(function(d) {
      dataset = d;
    }));

    $.when.apply($, promises).done(function () {
      var writer = new phantasus.GctWriter();

      var fileName = 'test.gct';

      var blobs = [];
      var textArray = [];
      var proxy = {
        push: function (text) {
          textArray.push(text);
          if (textArray.length === 10000) {
            var blob = new Blob([textArray.join('')], {type: 'text/plain;charset=charset=utf-8'});
            textArray = [];
            blobs.push(blob);
          }
        },
        join: function () {
          if (textArray.length > 0) {
            var blob = new Blob([textArray.join('')], {type: 'text/plain;charset=charset=utf-8'});
            blobs.push(blob);
            textArray = [];
          }

          var blob = new Blob(blobs, {type: 'text/plain;charset=charset=utf-8'});
          saveAs(blob, fileName, true);
        }
      };

      expect(writer.write(dataset, proxy)).not.toBeUndefined();

      done();
    });

  });
});