describe('dataset_test', function() {

	describe('transpose', function() {
		it('transposes a square matrix', function() {
			expect(new phantasus.TransposedDatasetView(new phantasus.Dataset({
				array : [ [ 1, 2 ], [ 3, 4 ] ],
				rows : 2,
				columns : 2
			}))).toBeDatasetValues(new phantasus.Dataset({
				array : [ [ 1, 3 ], [ 2, 4 ] ],
				rows : 2,
				columns : 2
			}));

		});

		it('transposes a non-square matrix', function() {
			expect(new phantasus.TransposedDatasetView(new phantasus.Dataset({
				array : [ [ 1, 2, 3, 4, 5 ], [ 2, 4, 6, 8, 10 ] ],
				rows : 2,
				columns : 5
			}))).toBeDatasetValues(new phantasus.Dataset({
				array : [ [ 1, 2 ], [ 2, 4 ], [ 3, 6 ], [ 4, 8 ], [ 5, 10 ] ],
				rows : 5,
				columns : 2
			}));

		});
		it('transposes a single-row matrix', function() {
			expect(new phantasus.TransposedDatasetView(new phantasus.Dataset({
				array : [ [ 1, 2, 3, 4, 5 ] ],
				rows : 1,
				columns : 5
			}))).toBeDatasetValues(new phantasus.Dataset({
				array : [ [ 1 ], [ 2 ], [ 3 ], [ 4 ], [ 5 ] ],
				rows : 5,
				columns : 1
			}));

		});

	});

	describe('slice', function() {
		it('slice rows', function() {
			expect(new phantasus.SlicedDatasetView(new phantasus.Dataset({
				array : [ [ 1, 2 ], [ 3, 4 ] ],
				rows : 2,
				columns : 2
			}), [ 0 ], null)).toBeDatasetValues(new phantasus.Dataset({
				array : [ [ 1, 2 ] ],
				rows : 1,
				columns : 2
			}));

		});

		it('slice columns', function() {
			expect(new phantasus.SlicedDatasetView(new phantasus.Dataset({
				array : [ [ 1, 2 ], [ 3, 4 ] ],
				rows : 2,
				columns : 2
			}), null, [ 0 ])).toBeDatasetValues(new phantasus.Dataset({
				array : [ [ 1 ], [ 3 ] ],
				rows : 2,
				columns : 1
			}));

		});

	});

});
