phantasus.HeatMapTrackShapeLegend = function (tracks, shapeModel) {
  phantasus.AbstractCanvas.call(this, false);
  this.tracks = tracks;
  this.shapeModel = shapeModel;
  this.canvas.style.position = '';
};
phantasus.HeatMapTrackShapeLegend.prototype = {
  getPreferredSize: function () {
    var tracks = this.tracks;
    var shapeModel = this.shapeModel;
    var canvas = this.canvas;
    var context = canvas.getContext('2d');
    var xpix = 0;
    var ypix = 0;
    var maxYPix = 0;
    for (var i = 0; i < tracks.length; i++) {
      ypix = 0;
      var maxWidth = 0;
      var vector = tracks[i].getVector();
      var map = shapeModel.getMap(vector.getName());

      map.forEach(function (color, key) {
        var width = context.measureText(key).width;
        if (!isNaN(width)) {
          maxWidth = Math.max(maxWidth, width);
        }
        ypix += 14;
      });

      xpix += maxWidth + 24;
      maxYPix = Math.max(maxYPix, ypix);
    }
    return {
      width: xpix,
      height: maxYPix > 0 ? (maxYPix + 30) : 0
    };
  },
  draw: function (clip, context) {
    // draw legends horizontally
    var tracks = this.tracks;
    var shapeModel = this.shapeModel;
    var xpix = 0;
    var ypix = 0;
    context.textAlign = 'left';
    context.textBaseline = 'top';
    context.font = '12px ' + phantasus.CanvasUtil.FONT_NAME;
    context.fillStyle = phantasus.CanvasUtil.FONT_COLOR;
    context.strokeStyle = 'black';
    for (var i = 0; i < tracks.length; i++) {
      ypix = 0;
      var maxWidth = 0;
      var vector = tracks[i].getVector();
      context.fillText(vector.getName(), xpix, ypix);
      maxWidth = Math.max(maxWidth,
        context.measureText(vector.getName()).width);
      ypix += 14;
      var map = shapeModel.getMap(vector.getName());
      var values = map.keys().sort(phantasus.SortKey.ASCENDING_COMPARATOR);
      values.forEach(function (key) {
        var shape = shapeModel.getMappedValue(vector, key);
        var width = context.measureText(key).width;
        if (!isNaN(width)) {
          maxWidth = Math.max(maxWidth, width);
        }
        phantasus.CanvasUtil.drawShape(context, shape, xpix + 8,
          ypix + 6, 6);
        context.fillText(key, xpix + 16, ypix);
        ypix += 14;
      });

      xpix += maxWidth + 24; // space between columns + shape
    }
  }
};
phantasus.Util.extend(phantasus.HeatMapTrackShapeLegend, phantasus.AbstractCanvas);
