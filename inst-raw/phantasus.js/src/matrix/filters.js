phantasus.CombinedFilter = function (isAndFilter) {
  this.filters = [];
  this.isAndFilter = isAndFilter;
  this.enabledFilters = [];
  this.name = 'combined filter';
};

phantasus.CombinedFilter.prototype = {
  shallowClone: function () {
    var f = new phantasus.CombinedFilter(this.isAndFilter);
    f.filters = this.filters.slice(0);
    return f;
  },
  isColumns: function () {
    return this.filters[0].isColumns();
  },
  toString: function () {
    return this.name;
  },
  setAnd: function (isAndFilter, notify) {
    this.isAndFilter = isAndFilter;
    if (notify) {
      this.trigger('and', {});
    }
  },
  isAnd: function () {
    return this.isAndFilter;
  },
  equals: function (f) {
    if (!(f instanceof phantasus.CombinedFilter)) {
      return false;
    }
    if (this.isAndFilter !== f.isAndFilter) {
      return false;
    }
    if (this.filters.length !== f.filters.length) {
      return false;
    }
    for (var i = 0, length = this.filters.length; i < length; i++) {
      if (!this.filters[i].equals(f.filters[i])) {
        return false;
      }
    }
    return true;
  },
  add: function (filter, notify) {
    this.filters.push(filter);
    if (notify) {
      this.trigger('add', {
        filter: filter
      });
    }
  },
  getFilters: function () {
    return this.filters;
  },
  get: function (index) {
    return this.filters[index];
  },
  indexOf: function (name, type) {
    for (var i = 0, length = this.filters.length; i < length; i++) {
      if (this.filters[i].toString() === name
        && (type == null ? true : this.filters[i] instanceof type)) {
        return i;
      }
    }
    return -1;
  },
  remove: function (index, notify) {
    this.filters.splice(index, 1);
    if (notify) {
      this.trigger('remove', {
        index: index
      });
    }
  },
  set: function (index, filter) {
    this.filters[index] = filter;
  },
  insert: function (index, filter) {
    this.filters.splice(index, 0, filter);
  },
  clear: function () {
    this.filters = [];
  },
  init: function (dataset) {
    for (var i = 0, nfilters = this.filters.length; i < nfilters; i++) {
      if (this.filters[i].isColumns()) { // all filters operate on rows
        this.filters[i].init(new phantasus.TransposedDatasetView(dataset));
      } else {
        this.filters[i].init(dataset);
      }

    }
    this.enabledFilters = this.filters.filter(function (filter) {
      return filter.isEnabled();
    });
  },
  accept: function (index) {
    var filters = this.enabledFilters;
    if (this.isAndFilter) {
      for (var i = 0, nfilters = filters.length; i < nfilters; i++) {
        if (filters[i].accept(index) === false) {
          return false;
        }
      }
      return true;
    } else {
      for (var i = 0, nfilters = filters.length; i < nfilters; i++) {
        if (filters[i].accept(index)) {
          return true;
        }
      }
      return false;
    }
  },
  isEnabled: function () {
    return this.enabledFilters.length > 0;
  }
};
phantasus.Util.extend(phantasus.CombinedFilter, phantasus.Events);
/**
 * @param acceptIndicesSet
 *            a phantasus.Set that contains the model indices in the dataset to
 *            retain.
 */
phantasus.IndexFilter = function (acceptIndicesSet, name, isColumns) {
  this.acceptIndicesSet = acceptIndicesSet;
  this.name = name;
  this.columns = isColumns;
};
phantasus.IndexFilter.prototype = {
  enabled: true,
  isColumns: function () {
    return this.columns;
  },
  isEnabled: function () {
    return this.enabled;
  },
  setAcceptIndicesSet: function (acceptIndicesSet) {
    this.acceptIndicesSet = acceptIndicesSet;
  },
  setEnabled: function (enabled) {
    this.enabled = enabled;
  },
  equals: function (filter) {
    return filter instanceof phantasus.IndexFilter
      && this.acceptIndicesSet.equals(filter.acceptIndicesSet);
  },
  init: function (dataset) {
  },
  toString: function () {
    return this.name;
  },
  /**
   *
   * @param index
   *            The model index in the dataset
   * @returns {Boolean} true if index passes filter
   */
  accept: function (index) {
    return this.acceptIndicesSet.has(index);
  }
};
phantasus.VectorFilter = function (set, maxSetSize, name, isColumns) {
  this.set = set;
  this.name = name;
  this.maxSetSize = maxSetSize;
  this.columns = isColumns;
};

phantasus.VectorFilter.prototype = {
  enabled: true,
  isColumns: function () {
    return this.columns;
  },
  isEnabled: function () {
    return this.enabled && this.set.size() > 0
      && this.set.size() !== this.maxSetSize && this.vector != null;
  },
  setEnabled: function (enabled) {
    this.enabled = enabled;
  },
  equals: function (filter) {
    return filter instanceof phantasus.VectorFilter
      && this.name === filter.name;
  },
  init: function (dataset) {
    this.vector = dataset.getRowMetadata().getByName(this.name);
  },
  toString: function () {
    return this.name;
  },
  /**
   *
   * @param index
   *            The model index in the dataset
   * @returns {Boolean} true if index passes filter
   */
  accept: function (index) {
    return this.set.has(this.vector.getValue(index));
  }
};

phantasus.NotNullFilter = function (name, isColumns) {
  this.name = name;
  this.columns = isColumns;
};
phantasus.NotNullFilter.prototype = {
  enabled: true,
  isColumns: function () {
    return this.columns;
  },
  isEnabled: function () {
    return this.enabled && this.vector != null;
  },
  setEnabled: function (enabled) {
    this.enabled = enabled;
  },
  equals: function (filter) {
    return filter instanceof phantasus.NotNullFilter
      && this.name === filter.name;
  },
  init: function (dataset) {
    this.vector = dataset.getRowMetadata().getByName(this.name);
  },
  toString: function () {
    return this.name;
  },
  /**
   *
   * @param index
   *            The model index in the dataset
   * @returns {Boolean} true if index passes filter
   */
  accept: function (index) {
    return this.vector.getValue(index) != null;
  }
};

phantasus.RangeFilter = function (min, max, name, isColumns) {
  this.min = min;
  this.max = max;
  this.name = name;
  this.columns = isColumns;
};

phantasus.RangeFilter.prototype = {
  enabled: true,
  isColumns: function () {
    return this.columns;
  },
  isEnabled: function () {
    return this.enabled && (!isNaN(this.min) || !isNaN(this.max))
      && this.vector;
  },
  setEnabled: function (enabled) {
    this.enabled = enabled;
  },
  setMin: function (value) {
    this.min = isNaN(value) ? -Number.MAX_VALUE : value;
  },
  setMax: function (value) {
    this.max = isNaN(value) ? Number.MAX_VALUE : value;
  },
  equals: function (filter) {
    return filter instanceof phantasus.RangeFilter
      && this.name === filter.name;
  },
  init: function (dataset) {
    this.vector = dataset.getRowMetadata().getByName(this.name);

  },
  toString: function () {
    return this.name;
  },
  /**
   *
   * @param index
   *            The model index in the dataset
   * @returns {Boolean} true if index passes filter
   */
  accept: function (index) {
    var value = this.vector.getValue(index);
    return value >= this.min && value <= this.max;
  }
};

phantasus.TopNFilter = function (n, direction, name, isColumns) {
  this.n = n;
  this.direction = direction;
  this.name = name;
  this.columns = isColumns;
};

phantasus.TopNFilter.TOP = 0;
phantasus.TopNFilter.BOTTOM = 1;
phantasus.TopNFilter.TOP_BOTTOM = 2;
phantasus.TopNFilter.prototype = {
  enabled: true,
  isColumns: function () {
    return this.columns;
  },
  isEnabled: function () {
    return this.enabled && this.n > 0 && this.vector;
  },
  setEnabled: function (enabled) {
    this.enabled = enabled;
  },
  setN: function (value) {
    this.n = value;
  },
  /**
   *
   * @param direction
   *            one of '
   */
  setDirection: function (direction) {
    this.direction = direction;
  },
  equals: function (filter) {
    return filter instanceof phantasus.TopNFilter
      && this.name === filter.name && this.n === filter.n
      && this.direction === filter.direction;
  },

  init: function (dataset) {
    if (!this.vector) {
      var vector = dataset.getRowMetadata().getByName(this.name);
      this.vector = vector;
      var set = new phantasus.Set();
      for (var i = 0, size = vector.size(); i < size; i++) {
        var value = vector.getValue(i);
        if (!isNaN(value)) {
          set.add(value);
        }
      }
      var values = set.values();
      // ascending order
      values.sort(function (a, b) {
        return (a === b ? 0 : (a < b ? -1 : 1));
      });
      this.sortedValues = values;
    }
    var topAndBottomIndices = [(this.sortedValues.length - this.n),
      (this.n - 1)];

    for (var i = 0; i < topAndBottomIndices.length; i++) {
      topAndBottomIndices[i] = Math.max(0, topAndBottomIndices[i]);
      topAndBottomIndices[i] = Math.min(this.sortedValues.length - 1,
        topAndBottomIndices[i]);
    }

    var topAndBottomValues = [this.sortedValues[topAndBottomIndices[0]],
      this.sortedValues[topAndBottomIndices[1]]];

    if (this.direction === phantasus.TopNFilter.TOP) {
      this.f = function (val) {
        return isNaN(val) ? false : val >= topAndBottomValues[0];
      };
    } else if (this.direction === phantasus.TopNFilter.BOTTOM) {
      this.f = function (val) {
        return isNaN(val) ? false : val <= topAndBottomValues[1];
      };
    } else {
      this.f = function (val) {
        return isNaN(val) ? false
          : (val >= topAndBottomValues[0] || val <= topAndBottomValues[1]);
      };
    }

  },
  /**
   *
   * @param index
   *            The model index in the dataset
   * @returns {Boolean} true if index passes filter
   */
  accept: function (index) {
    return this.f(this.vector.getValue(index));
  },
  toString: function () {
    return this.name;
  }
};

phantasus.AlwaysTrueFilter = function () {

};

phantasus.AlwaysTrueFilter.prototype = {
  isEnabled: function () {
    return false;
  },
  setEnabled: function (enabled) {

  },
  equals: function (filter) {
    return filter instanceof phantasus.AlwaysTrueFilter;

  },
  init: function (dataset) {

  },
  toString: function () {
    return 'AlwaysTrue';
  },
  /**
   *
   * @param index
   *            The model index in the dataset
   * @returns {Boolean} true if index passes filter
   */
  accept: function (index) {
    return true;
  }
};

phantasus.CombinedFilter.fromJSON = function (combinedFilter, json) {
  combinedFilter.setAnd(json.isAnd);
  json.filters.forEach(function (filter) {
    if (filter.type === 'set') {
      var set = new phantasus.Set();
      filter.values.forEach(function (value) {
        set.add(value);
      });
      combinedFilter.add(new phantasus.VectorFilter(
        set,
        filter.maxSetSize,
        filter.name,
        filter.isColumns
      ));
    } else if (filter.type === 'range') {
      combinedFilter.add(new phantasus.RangeFilter(
        filter.min,
        filter.max,
        filter.name,
        filter.isColumns
      ));
    } else if (filter.type === 'top') {
      combinedFilter.add(new phantasus.TopNFilter(
        filter.n,
        filter.direction,
        filter.name,
        filter.isColumns
      ));
    } else if (filter.type === 'index') {
      var set = new phantasus.Set();
      filter.indices.forEach(function (value) {
        set.add(value);
      });
      combinedFilter.add(new phantasus.IndexFilter(
        set,
        filter.name,
        filter.isColumns
      ));
    } else {
      // console.log('Unknown filter type');
    }
  });
};

phantasus.CombinedFilter.toJSON = function (filter) {
  var json = {
    isAnd: filter.isAnd(),
    filters: []
  };
  filter.getFilters().forEach(function (filter) {
    if (filter.isEnabled()) {
      if (filter instanceof phantasus.VectorFilter) {
        // phantasus.VectorFilter = function (set, maxSetSize, name, isColumns)
        json.filters.push({
          name: filter.name,
          isColumns: filter.isColumns(),
          values: filter.set.values(),
          maxSetSize: filter.maxSetSize,
          type: 'set'
        });
      } else if (filter instanceof phantasus.RangeFilter) {
        // phantasus.RangeFilter = function (min, max, name, isColumns)
        json.filters.push({
          name: filter.name,
          isColumns: filter.isColumns(),
          min: filter.min,
          max: filter.max,
          type: 'range'
        });
      } else if (filter instanceof phantasus.TopNFilter) {
        // phantasus.TopNFilter = function (n, direction, name, isColumns)
        json.filters.push({
          name: filter.name,
          isColumns: filter.isColumns(),
          min: filter.n,
          max: filter.direction,
          type: 'top'
        });
      } else if (filter instanceof phantasus.IndexFilter) {
        // phantasus.IndexFilter = function (acceptIndicesSet, name, isColumns
        json.filters.push({
          name: filter.name,
          isColumns: filter.isColumns(),
          indices: filter.acceptIndicesSet.values(),
          type: 'index'
        });
      }
    }
  });
  return json;
};
