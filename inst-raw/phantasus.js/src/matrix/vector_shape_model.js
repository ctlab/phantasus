phantasus.VectorShapeModel = function () {
  this.shapes = phantasus.VectorShapeModel.SHAPES;
  this.vectorNameToShapeMap = new phantasus.Map();
};

phantasus.VectorShapeModel.SHAPES = ['circle', 'square', 'plus', 'x',
  'asterisk', 'diamond', 'triangle-up', 'triangle-down', 'triangle-left',
  'triangle-right', 'minus'];
phantasus.VectorShapeModel.STANDARD_SHAPES = {
  cp: 'diamond',
  oe: 'plus',
  pcl: 'asterisk',
  kd: 'minus',
  ctrl: 'circle'
};

phantasus.VectorShapeModel.prototype = {
  clear: function (vector) {
    this.vectorNameToShapeMap.remove(vector.getName());
  },
  copy: function () {
    var c = new phantasus.VectorShapeModel();
    c.shapes = this.shapes.slice(0);
    this.vectorNameToShapeMap.forEach(function (shapeMap, name) {
      var newShapeMap = new phantasus.Map();
      newShapeMap.setAll(shapeMap); // copy existing values
      c.vectorNameToShapeMap.set(name, newShapeMap);
    });

    return c;
  },
  clearAll: function () {
    this.vectorNameToShapeMap = new phantasus.Map();
  },
  _getShapeForValue: function (value) {
    if (value == null) {
      return 'none';
    }
    var str = value.toString().toLowerCase();
    var mapped = phantasus.VectorShapeModel.STANDARD_SHAPES[str];
    if (mapped !== undefined) {
      return mapped;
    }

    // try to reuse existing map
    var existingMetadataValueToShapeMap = this.vectorNameToShapeMap
      .values();
    for (var i = 0, length = existingMetadataValueToShapeMap.length; i < length; i++) {
      var shape = existingMetadataValueToShapeMap[i].get(value);
      if (shape !== undefined) {
        return shape;
      }
    }

  },
  getMap: function (name) {
    return this.vectorNameToShapeMap.get(name);
  },
  getMappedValue: function (vector, value) {
    var metadataValueToShapeMap = this.vectorNameToShapeMap.get(vector
      .getName());
    if (metadataValueToShapeMap === undefined) {
      metadataValueToShapeMap = new phantasus.Map();
      this.vectorNameToShapeMap.set(vector.getName(),
        metadataValueToShapeMap);
      // set all possible shapes
      var values = phantasus.VectorUtil.getValues(vector);
      for (var i = 0, nvalues = values.length; i < nvalues; i++) {
        var shape = this._getShapeForValue(values[i]);
        if (shape == null) {
          shape = this.shapes[i % this.shapes.length];
        }
        metadataValueToShapeMap.set(values[i], shape);
      }
    }
    var shape = metadataValueToShapeMap.get(value);
    if (shape == null) {
      shape = this._getShapeForValue(value);
      if (shape == null) {
        var index = metadataValueToShapeMap.size();
        shape = this.shapes[index % this.shapes.length];
      }
      metadataValueToShapeMap.set(value, shape);
    }
    return shape;
  },
  setMappedValue: function (vector, value, shape) {
    var metadataValueToShapeMap = this.vectorNameToShapeMap.get(vector
      .getName());
    if (metadataValueToShapeMap === undefined) {
      metadataValueToShapeMap = new phantasus.Map();
      this.vectorNameToShapeMap.set(vector.getName(),
        metadataValueToShapeMap);
    }
    metadataValueToShapeMap.set(value, shape);
  }
};
