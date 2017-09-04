phantasus.DiscreteColorSupplier = function () {
  this.colorMap = new phantasus.Map();
  this.hiddenValue = 0;
  this.hiddenValues = new phantasus.Set();
  phantasus.AbstractColorSupplier.call(this);
  this.scalingMode = phantasus.HeatMapColorScheme.ScalingMode.FIXED;
};

phantasus.DiscreteColorSupplier.prototype = {
  createInstance: function () {
    return new phantasus.DiscreteColorSupplier();
  },
  /**
   * @param.array Array of name, value, color pairs
   */
  setColorMap: function (array) {
    this.colorMap = new phantasus.Map();
    this.colors = [];
    this.fractions = [];
    this.names = [];
    this.min = Number.MAX_VALUE;
    this.max = -Number.MAX_VALUE;
    for (var i = 0; i < array.length; i++) {
      this.colorMap.set(array[i].value, array[i].color);
      this.fractions.push(array[i].value);
      this.names.push(array[i].name);
      this.colors.push(array[i].color);
      this.min = Math.min(this.min, array[i].value);
      this.max = Math.max(this.max, array[i].value);
    }
  },
  copy: function () {
    var c = this.createInstance();
    c.names = this.names.slice(0);
    c.colorMap = new phantasus.Map();
    this.colorMap.forEach(function (color, value) {
      c.colorMap.set(value, color);
    });
    c.colors = this.colors.slice(0);
    c.fractions = this.fractions.slice(0);
    this.hiddenValues.forEach(function (val) {
      c.hiddenValues.add(val);
    });

    c.missingColor = this.missingColor;
    return c;
  },

  isStepped: function () {
    return true;
  },
  getColor: function (row, column, value) {
    if (this.hiddenValues.has(value)) {
      value = this.hiddenValue;
    }

    if (isNaN(value)) {
      return this.missingColor;
    }
    return this.colorMap.get(value);
  }
};
phantasus.Util.extend(phantasus.DiscreteColorSupplier,
  phantasus.AbstractColorSupplier);
