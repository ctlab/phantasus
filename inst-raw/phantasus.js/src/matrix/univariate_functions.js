/**
 * Provides percentile computation.
 * <p>
 * There are several commonly used methods for estimating percentiles (a.k.a.
 * quantiles) based on sample data. For large samples, the different methods
 * agree closely, but when sample sizes are small, different methods will give
 * significantly different results. The algorithm implemented here works as
 * follows:
 * <ol>
 * <li>Let <code>n</code> be the length of the (sorted) array and
 * <code>0 < p <= 100</code> be the desired percentile.</li>
 * <li>If <code> n = 1 </code> return the unique array element (regardless of
 * the value of <code>p</code>); otherwise</li>
 * <li>Compute the estimated percentile position
 * <code> pos = p * (n + 1) / 100</code> and the difference, <code>d</code>
 * between <code>pos</code> and <code>floor(pos)</code> (i.e. the fractional
 * part of <code>pos</code>). If <code>pos >= n</code> return the largest
 * element in the array; otherwise</li>
 * <li>Let <code>lower</code> be the element in position
 * <code>floor(pos)</code> in the array and let <code>upper</code> be the
 * next element in the array. Return <code>lower + d * (upper - lower)</code></li>
 * </ol>
 *
 * @param p Percentile between 0 and 100
 */
phantasus.Percentile = function (vector, p, isSorted) {
  return phantasus.ArrayPercentile(phantasus.RemoveNaN(vector), p, isSorted);
};
/**
 * @private
 * @ignore
 */
phantasus.RemoveNaN = function (values) {
  var array = [];
  for (var i = 0, size = values.size(); i < size; i++) {
    var value = values.getValue(i);
    if (!isNaN(value)) {
      array.push(value);
    }
  }
  return array;
};
phantasus.Median = function (vector) {
  return phantasus.ArrayPercentile(phantasus.RemoveNaN(vector), 50, false);
};
phantasus.Median.toString = function () {
  return 'Median';
};
/**
 * @ignore
 */
phantasus.ArrayPercentile = function (values, p, isSorted) {

  if (!isSorted) {
    values.sort(function (a, b) {
      return (a < b ? -1 : (a === b ? 0 : 1));
    });
  }
  return d3.quantile(values, p / 100);
};
/**
 * @ignore
 */
phantasus.MaxPercentiles = function (percentiles) {
  var f = function (vector) {
    var values = [];
    for (var i = 0, size = vector.size(); i < size; i++) {
      var value = vector.getValue(i);
      if (!isNaN(value)) {
        values.push(value);
      }
    }
    if (values.length === 0) {
      return NaN;
    }
    values.sort(function (a, b) {
      return (a < b ? -1 : (a === b ? 0 : 1));
    });
    var max = 0;
    for (var i = 0; i < percentiles.length; i++) {
      var p = phantasus.ArrayPercentile(values, percentiles[i], true);
      if (Math.abs(p) > Math.abs(max)) {
        max = p;
      }
    }
    return max;
  };
  f.toString = function () {
    var s = ['Maximum of '];
    for (var i = 0, length = percentiles.length; i < length; i++) {
      if (i > 0 && length > 2) {
        s.push(', ');
      }
      if (i === length - 1) {
        s.push(length == 2 ? ' and ' : 'and ');
      }
      s.push(percentiles[i]);
    }
    s.push(' percentiles');
    return s.join('');
  };
  return f;
};
phantasus.Mean = function (vector) {
  var sum = 0;
  var count = 0;
  for (var i = 0, length = vector.size(); i < length; i++) {
    var val = vector.getValue(i);
    if (!isNaN(val)) {
      sum += val;
      count++;
    }
  }
  return count === 0 ? NaN : sum / count;
};
phantasus.Mean.toString = function () {
  return 'Mean';
};
phantasus.Sum = function (vector) {
  var sum = 0;
  var found = false;
  for (var i = 0, length = vector.size(); i < length; i++) {
    var val = vector.getValue(i);
    if (!isNaN(val)) {
      found = true;
      sum += val;
    }
  }
  return !found ? NaN : sum;
};
phantasus.Sum.toString = function () {
  return 'Sum';
};
phantasus.CountNonNaN = function (vector) {
  var count = 0;
  for (var i = 0, length = vector.size(); i < length; i++) {
    var val = vector.getValue(i);
    if (!isNaN(val)) {
      count++;
    }
  }
  return count;
};
phantasus.CountNonNaN.toString = function () {
  return 'Count non-NaN';
};

phantasus.Max = function (vector) {
  var max = -Number.MAX_VALUE;
  var found = false;
  for (var i = 0, length = vector.size(); i < length; i++) {
    var val = vector.getValue(i);
    if (!isNaN(val)) {
      found = true;
      max = Math.max(max, val);
    }
  }
  return !found ? NaN : max;
};
phantasus.Max.toString = function () {
  return 'Max';
};
phantasus.Min = function (vector) {
  var min = Number.MAX_VALUE;
  var found = false;
  for (var i = 0, length = vector.size(); i < length; i++) {
    var val = vector.getValue(i);
    if (!isNaN(val)) {
      found = true;
      min = Math.min(min, val);
    }
  }
  return !found ? NaN : min;
};
phantasus.Min.toString = function () {
  return 'Min';
};
phantasus.Variance = function (list, mean) {
  if (mean == undefined) {
    mean = phantasus.Mean(list);
  }
  var sum = 0;
  var n = 0;
  for (var j = 0, size = list.size(); j < size; j++) {
    var x = list.getValue(j);
    if (!isNaN(x)) {
      var diff = x - mean;
      diff = diff * diff;
      sum += diff;
      n++;
    }
  }
  if (n <= 1) {
    return NaN;
  }
  n = n - 1;
  if (n < 1) {
    n = 1;
  }
  var variance = sum / n;
  return variance;
};
phantasus.Variance.toString = function () {
  return 'Variance';
};

phantasus.StandardDeviation = function (list, mean) {
  return Math.sqrt(phantasus.Variance(list, mean));
};
phantasus.StandardDeviation.toString = function () {
  return 'Standard deviation';
};

var LOG_10 = Math.log(10);
phantasus.Log10 = function (x) {
  return x <= 0 ? 0 : Math.log(x) / LOG_10;
};
var LOG_2 = Math.log(2);
phantasus.Log2 = function (x) {
  return x <= 0 ? 0 : Math.log(x) / LOG_2;
};

/**
 * Computes the False Discovery Rate using the BH procedure.
 *
 * @param nominalPValues
 *            Array of nominal p-values.
 */
phantasus.FDR_BH = function (nominalPValues) {
  var size = nominalPValues.length;
  var fdr = [];
  var pValueIndices = phantasus.Util.indexSort(nominalPValues, true);
  var ranks = phantasus.Util.rankIndexArray(pValueIndices);

  // check for ties
  for (var i = pValueIndices.length - 1; i > 0; i--) {
    var bigPValue = nominalPValues[pValueIndices[i]];
    var smallPValue = nominalPValues[pValueIndices[i - 1]];
    if (bigPValue == smallPValue) {
      ranks[pValueIndices[i - 1]] = ranks[pValueIndices[i]];
    }
  }
  for (var i = 0; i < size; i++) {
    var rank = ranks[i];
    var p = nominalPValues[i];
    fdr[i] = (p * size) / rank;
  }

  // ensure fdr is monotonically decreasing
  var pIndices = phantasus.Util.indexSort(nominalPValues, false);
  for (var i = 0; i < pIndices.length - 1; i++) {
    var highIndex = pIndices[i];
    var lowIndex = pIndices[i + 1];
    fdr[lowIndex] = Math.min(fdr[lowIndex], fdr[highIndex]);
  }
  for (var i = 0; i < size; i++) {
    fdr[i] = Math.min(fdr[i], 1);
  }
  return fdr;
};

phantasus.FDR_BH.tString = function () {
  return 'FDR(BH)';
};

phantasus.MAD = function (list, median) {
  if (median == null) {
    median = phantasus.Percentile(list, 50);
  }
  var temp = [];
  for (var j = 0, size = list.size(); j < size; j++) {
    var value = list.getValue(j);
    if (!isNaN(value)) {
      temp.push(Math.abs(value - median));
    }
  }
  var r = phantasus.Percentile(new phantasus.Vector('', temp.length)
    .setArray(temp), 50);
  return 1.4826 * r;
};
phantasus.MAD.toString = function () {
  return 'Median absolute deviation';
};
phantasus.CV = function (list) {
  var mean = phantasus.Mean(list);
  var stdev = Math.sqrt(phantasus.Variance(list, mean));
  return stdev / mean;
};
phantasus.CV.toString = function () {
  return 'Coefficient of variation';
};

phantasus.BoxPlotItem = function (list) {
  var values = phantasus.RemoveNaN(list);
  values.sort(function (a, b) {
    return (a === b ? 0 : (a < b ? -1 : 1));
  });
  if (values.length === 0) {
    return {
      median: NaN,
      q1: NaN,
      q3: NaN,
      lowerAdjacentValue: NaN,
      upperAdjacentValue: NaN
    };
  }
  var median = phantasus.ArrayPercentile(values, 50, true);
  var q1 = phantasus.ArrayPercentile(values, 25, true);
  var q3 = phantasus.ArrayPercentile(values, 75, true);
  var w = 1.5;
  var upperAdjacentValue = -Number.MAX_VALUE;
  var lowerAdjacentValue = Number.MAX_VALUE;
  // The upper adjacent value (UAV) is the largest observation that is
  // less than or equal to
  // the upper inner fence (UIF), which is the third quartile plus
  // 1.5*IQR.
  //
  // The lower adjacent value (LAV) is the smallest observation that is
  // greater than or equal
  // to the lower inner fence (LIF), which is the first quartile minus
  // 1.5*IQR.
  var upperOutlier = q3 + w * (q3 - q1);
  var lowerOutlier = q1 - w * (q3 - q1);
  var sum = 0;
  for (var i = 0, length = values.length; i < length; i++) {
    var value = values[i];
    if (value <= upperOutlier) {
      upperAdjacentValue = Math.max(upperAdjacentValue, value);
    }
    if (value >= lowerOutlier) {
      lowerAdjacentValue = Math.min(lowerAdjacentValue, value);
    }
    sum += value;
    // if (value > upperOutlier) {
    // upperOutliers.add(new Outlier(i, j, value));
    // }
    // if (value < lowerOutlier) {
    // lowerOutliers.add(new Outlier(i, j, value));
    // }
  }
  var mean = sum / values.length;
  if (lowerAdjacentValue > q1) {
    lowerAdjacentValue = q1;
  }
  if (upperAdjacentValue < q3) {
    upperAdjacentValue = q3;
  }

  return {
    mean: mean,
    median: median,
    q1: q1, // Lower Quartile
    q3: q3, // Upper Quartile
    lowerAdjacentValue: lowerAdjacentValue, // Lower Whisker
    upperAdjacentValue: upperAdjacentValue
    // Upper Whisker
  };

};
