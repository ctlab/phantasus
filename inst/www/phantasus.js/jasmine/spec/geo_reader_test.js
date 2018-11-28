jasmine.DEFAULT_TIMEOUT_INTERVAL = 240000;

describe('geo_reader_test', function () {

  describe('loading mono gse datasets', function () {
    var result;
    beforeEach(function (done) {
      var reader = new phantasus.GeoReader();
      reader.read('GSE53986', function (err, success) {
        result = success;

        done();
      })
    });

    it('loades dataset GSE53986 correctly', function () {
      expect(result).not.toBeUndefined();
      expect(result.length).toEqual(1);

      var dataset = result[0];

      expect(dataset.getRowCount() === 45101 || dataset.getRowCount() === 45100).toEqual(true);
      expect(dataset.getColumnCount()).toEqual(16);
    });
  });

});
