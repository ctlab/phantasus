/**
 * Creates a new meta data model instance.
 *
 * @param itemCount {number}
 *            the number of items that vectors in this instances will hold.
 * @implements {phantasus.MetadataModelInterface}
 * @constructor
 */
phantasus.MetadataModel = function (itemCount) {
  this.itemCount = itemCount;
  this.vectors = [];
};
phantasus.MetadataModel.prototype = {
  add: function (name, options) {
    var index = phantasus.MetadataUtil.indexOf(this, name);
    var oldVector;
    if (index !== -1) {
      oldVector = this.remove(index);
    }
    var v = new phantasus.Vector(name, this.getItemCount());
    if (oldVector != null) {
      // copy properties?
//			oldVector.getProperties().forEach(function(val, key) {
//				if (!phantasus.VectorKeys.COPY_IGNORE.has(key)) {
//					v.getProperties().set(key, val);
//				}
//
//			});
      // copy values
      for (var i = 0, size = oldVector.size(); i < size; i++) {
        var val = oldVector.getValue(i);
        v.setValue(i, val);
      }
    }
    this.vectors.push(v);
    return v;
  },
  getItemCount: function () {
    return this.itemCount;
  },
  get: function (index) {
    if (index < 0 || index >= this.vectors.length) {
      throw 'index ' + index + ' out of range';
    }
    return this.vectors[index];
  },
  remove: function (index) {
    if (index < 0 || index >= this.vectors.length) {
      throw 'index ' + index + ' out of range';
    }
    return this.vectors.splice(index, 1)[0];
  },
  getByName: function (name) {
    var index = phantasus.MetadataUtil.indexOf(this, name);
    return index !== -1 ? this.get(index) : undefined;
  },
  getMetadataCount: function () {
    return this.vectors.length;
  }
};
