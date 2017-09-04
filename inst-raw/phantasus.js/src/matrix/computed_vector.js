/**
 * Creates a new computed vector with the given name and size.
 *
 * @param name
 *            the vector name
 * @param size
 *            the number of elements in this vector
 * @param callback {Function} that takes an index and returns the value at the specified index
 * @constructor
 */
phantasus.ComputedVector = function (name, size, callback) {
  phantasus.AbstractVector.call(this, name, size);
  this.callback = callback;
};

phantasus.ComputedVector.prototype = {
  getValue: function (index) {
    return this.callback(index);
  }
};
phantasus.Util.extend(phantasus.ComputedVector, phantasus.AbstractVector);
