phantasus.ElementSelectionModel = function (project) {
  this.viewIndices = new phantasus.Set();
  this.project = project;
};
phantasus.ElementSelectionModel.prototype = {
  click: function (rowIndex, columnIndex, add) {
    var id = new phantasus.Identifier([rowIndex, columnIndex]);
    var isSelected = this.viewIndices.has(id);
    if (add) {
      isSelected ? this.viewIndices.remove(id) : this.viewIndices.add(id);
    } else {
      this.viewIndices.clear();
      if (!isSelected) {
        this.viewIndices.add(id);
      }
    }
    this.trigger('selectionChanged');
  },
  getProject: function () {
    return this.project;
  },
  setViewIndices: function (indices) {
    this.viewIndices = indices;
    this.trigger('selectionChanged');
  },
  clear: function () {
    this.viewIndices = new phantasus.Set();
  },
  /**
   *
   * @returns {phantasus.Set}
   */
  getViewIndices: function () {
    return this.viewIndices;
  },
  count: function () {
    return this.viewIndices.size();
  },
  toModelIndices: function () {
    var project = this.project;
    var modelIndices = [];
    this.viewIndices.forEach(function (id) {
      modelIndices.push(project
        .convertViewRowIndexToModel(id.getArray()[0]), project
        .convertViewColumnIndexToModel(id.getArray()[1]));
    });
    return modelIndices;
  },
  save: function () {
    this.modelIndices = this.toModelIndices();
  },
  restore: function () {
    var project = this.project;
    this.viewIndices = new phantasus.Set();
    for (var i = 0, length = this.modelIndices.length; i < length; i++) {
      var rowIndex = project
        .convertModelRowIndexToView(this.modelIndices[i][0]);
      var columnIndex = project
        .convertModelColumnIndexToView(this.modelIndices[i][1]);
      if (rowIndex !== -1 && columnIndex !== -1) {
        this.viewIndices.add(new phantasus.Identifier([rowIndex,
          columnIndex]));
      }
    }
  }
};
phantasus.Util.extend(phantasus.ElementSelectionModel, phantasus.Events);
