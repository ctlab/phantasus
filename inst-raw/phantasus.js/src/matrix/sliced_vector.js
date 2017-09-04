phantasus.SlicedVector = function (v, indices) {
  phantasus.VectorAdapter.call(this, v);
  this.indices = indices;
};
phantasus.SlicedVector.prototype = {
  setValue: function (i, value) {
    this.v.setValue(this.indices[i], value);
  },
  getValue: function (i) {
    return this.v.getValue(this.indices[i]);
  },
  size: function () {
    return this.indices.length;
  }
};
phantasus.Util.extend(phantasus.SlicedVector, phantasus.VectorAdapter);
