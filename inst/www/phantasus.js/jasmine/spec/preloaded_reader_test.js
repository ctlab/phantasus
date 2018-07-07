describe('preloaded_reader_test', function () {

  describe('loading mono preloaded datasets', function () {
    var result;

    beforeEach(function (done) {
      originalTimeout = jasmine.DEFAULT_TIMEOUT_INTERVAL;
      jasmine.DEFAULT_TIMEOUT_INTERVAL = 100000;

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

    afterEach(function() {
      jasmine.DEFAULT_TIMEOUT_INTERVAL = originalTimeout;
    });
  })
});
