/**
 * An ordered collection of values.
 *
 * Creates a new vector with the given name and size.
 *
 * @param name
 *            the vector name
 * @param size
 *            the number of elements in this vector
 * @constructor
 */
phantasus.Vector = function (name, size) {
  this.array = [];
  phantasus.AbstractVector.call(this, name, size);
};
/**
 * @static
 */
phantasus.Vector.fromArray = function (name, array) {
  var v = new phantasus.Vector(name, array.length);
  v.array = array;
  return v;
};
phantasus.Vector.prototype = {
  /**
   * @ignore
   * @param value
   */
  push: function (value) {
    this.array.push(value);
  },
  /**
   * Sets the value at the specified index.
   *
   * @param index
   *            the index
   * @param value
   *            the value
   */
  setValue: function (index, value) {
    this.array[index] = value;
  },
  getValue: function (index) {
    return this.array[index];
  },
  /**
   * @ignore
   * @param name
   */
  setName: function (name) {
    this.name = name;
  },
  /**
   * @ignore
   * @param array
   * @returns {phantasus.Vector}
   */
  setArray: function (array) {
    this.array = array;
    return this;
  }
};
phantasus.Util.extend(phantasus.Vector, phantasus.AbstractVector);
