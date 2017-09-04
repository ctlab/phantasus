var flag = false;

function testAsyncWithDeferred() {
  // Get a jQuery deferred
  var deferred = $.Deferred();

  // Wait two seconds, then set the flag to true
  setTimeout(function () {
    flag = true;

    // Resolve the deferred
    deferred.resolve();
  }, 2000);

  // Return the deferred promise
  return deferred.promise();
}

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

      expect(dataset.getRowCount()).toEqual(45101);
      expect(dataset.getColumnCount()).toEqual(16);
    });
  });

});