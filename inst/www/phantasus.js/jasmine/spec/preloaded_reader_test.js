describe('preloaded_reader_test', function () {

  describe('loading mono preloaded datasets', function () {
    var result;

    beforeEach(function (done) {
      var reader = new phantasus.PreloadedReader();
      reader.read('GSE27112-GPL6103', function (err, success) {
        result = success;

        done();
      })
    });

    it('loads dataset GSE27112-GPL6103 correctly', function () {
      expect(result).not.toBeUndefined();
      expect(result.length).toEqual(1);

      var dataset = result[0];

      expect(dataset.getRowCount()).toEqual(20);
      expect(dataset.getColumnCount()).toEqual(5);
    })
  })
});
