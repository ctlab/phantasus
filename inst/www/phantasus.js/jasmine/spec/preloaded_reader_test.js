describe('preloaded_reader_test', function () {

  describe('loading mono preloaded datasets', function () {
    var result;

    beforeEach(function (done) {
      var reader = new phantasus.PreloadedReader();
      reader.read('GSE53986', function (err, success) {
        result = success;

        done();
      })
    });

    it('loades dataset GSE53986 correctly', function () {
      expect(result).not.toBeUndefined();
      expect(result.length).toEqual(1);

      var dataset = result[0];

      expect(dataset.getRowCount()).toEqual(45101);
      expect(dataset.getColumnCount()).toEqual(16);
    })
  })
});