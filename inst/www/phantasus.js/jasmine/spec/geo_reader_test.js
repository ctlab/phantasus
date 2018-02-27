describe('geo_reader_test', function () {

  describe('loading mono gse datasets', function () {
    var originalTimeout;
    var result;
    beforeEach(function () {
      //Will run into default timeout on slow connections.
      originalTimeout = jasmine.DEFAULT_TIMEOUT_INTERVAL;
      jasmine.DEFAULT_TIMEOUT_INTERVAL = 600000;

    });

    it('loades dataset GSE53986 correctly', function (done) {
      var reader = new phantasus.GeoReader();
      reader.read('GSE53986', function (err, success) {
        result = success;
        expect(result).not.toBeUndefined();
        expect(result.length).toEqual(1);

        var dataset = result[0];

        expect(dataset.getRowCount()).toEqual(45101);
        expect(dataset.getColumnCount()).toEqual(16);

        done();
      })

    });

    afterEach(function() {
      jasmine.DEFAULT_TIMEOUT_INTERVAL = originalTimeout;
    });
  });

});
