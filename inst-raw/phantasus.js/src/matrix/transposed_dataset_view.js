phantasus.TransposedDatasetView = function (dataset) {
  phantasus.DatasetAdapter.call(this, dataset);
};
phantasus.TransposedDatasetView.prototype = {
  getRowCount: function () {
    return this.dataset.getColumnCount();
  },
  getColumnCount: function () {
    return this.dataset.getRowCount();
  },
  getValue: function (i, j, seriesIndex) {
    return this.dataset.getValue(j, i, seriesIndex);
  },
  setValue: function (i, j, value, seriesIndex) {
    this.dataset.setValue(j, i, value, seriesIndex);
  },
  getRowMetadata: function () {
    return this.dataset.getColumnMetadata();
  },
  getColumnMetadata: function () {
    return this.dataset.getRowMetadata();
  }
};
phantasus.Util.extend(phantasus.TransposedDatasetView, phantasus.DatasetAdapter);
