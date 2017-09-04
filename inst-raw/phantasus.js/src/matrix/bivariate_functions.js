phantasus.SignalToNoise = function (list1, list2) {
  var m1 = phantasus.Mean(list1);
  var m2 = phantasus.Mean(list2);
  var s1 = Math.sqrt(phantasus.Variance(list1, m1));
  var s2 = Math.sqrt(phantasus.Variance(list2, m2));
  return (m1 - m2) / (s1 + s2);
};
phantasus.SignalToNoise.toString = function () {
  return 'Signal to noise';
};

phantasus.createSignalToNoiseAdjust = function (percent) {
  percent = percent || 0.2;
  var f = function (list1, list2) {
    var m1 = phantasus.Mean(list1);
    var m2 = phantasus.Mean(list2);
    var s1 = Math.sqrt(phantasus.Variance(list1, m1));
    var s2 = Math.sqrt(phantasus.Variance(list2, m2));
    s1 = phantasus.SignalToNoise.thresholdStandardDeviation(m1, s1, percent);
    s2 = phantasus.SignalToNoise.thresholdStandardDeviation(m2, s2, percent);
    // ensure variance is at least 20% of mean
    return (m1 - m2) / (s1 + s2);
  };
  f.toString = function () {
    return 'Signal to noise (adjust standard deviation)';
  };
  return f;
};

phantasus.SignalToNoise.thresholdStandardDeviation = function (mean,
                                                              standardDeviation, percent) {
  var returnValue = standardDeviation;
  var absMean = Math.abs(mean);
  var minStdev = percent * absMean;
  if (minStdev > standardDeviation) {
    returnValue = minStdev;
  }

  if (returnValue < percent) {
    returnValue = percent;
  }
  return returnValue;
};

phantasus.createContingencyTable = function (listOne, listTwo, groupingValue) {
  if (groupingValue == null || isNaN(groupingValue)) {
    groupingValue = 1;
  }
  var aHit = 0;
  var aMiss = 0;
  for (var j = 0, size = listOne.size(); j < size; j++) {
    var val = listOne.getValue(j);
    if (!isNaN(val)) {
      if (val >= groupingValue) {
        aHit++;
      } else {
        aMiss++;
      }
    }

  }
  var bHit = 0;
  var bMiss = 0;
  for (var j = 0, size = listTwo.size(); j < size; j++) {
    var val = listTwo.getValue(j);
    if (!isNaN(val)) {
      if (val >= groupingValue) {
        bHit++;
      } else {
        bMiss++;
      }
    }

  }
  // listOne=drawn, listTwo=not drawn
  // green=1, red=0
  var N = aHit + aMiss + bHit + bMiss;
  var K = aHit + bHit;
  var n = aHit + aMiss;
  var k = aHit;
  var a = k;
  var b = K - k;
  var c = n - k;
  var d = N + k - n - K;
  return [a, b, c, d];
};
phantasus.FisherExact = function (listOne, listTwo) {
  var abcd = phantasus.createContingencyTable(listOne, listTwo, 1);
  return phantasus.FisherExact.fisherTest(abcd[0], abcd[1], abcd[2], abcd[3]);
};

phantasus.createFisherExact = function (groupingValue) {
  var f = function (listOne, listTwo) {
    var abcd = phantasus.createContingencyTable(listOne, listTwo,
      groupingValue);
    return phantasus.FisherExact.fisherTest(abcd[0], abcd[1], abcd[2],
      abcd[3]);
  };
  return f;

};

/**
 * Computes the hypergeometric probability.
 */
phantasus.FisherExact.phyper = function (a, b, c, d) {
  return Math
    .exp((phantasus.FisherExact.logFactorial(a + b)
      + phantasus.FisherExact.logFactorial(c + d)
      + phantasus.FisherExact.logFactorial(a + c) + phantasus.FisherExact
        .logFactorial(b + d))
      - (phantasus.FisherExact.logFactorial(a)
      + phantasus.FisherExact.logFactorial(b)
      + phantasus.FisherExact.logFactorial(c)
      + phantasus.FisherExact.logFactorial(d) + phantasus.FisherExact
        .logFactorial(a + b + c + d)));

};

phantasus.FisherExact.logFactorials = [0.00000000000000000,
  0.00000000000000000, 0.69314718055994531, 1.79175946922805500,
  3.17805383034794562, 4.78749174278204599, 6.57925121201010100,
  8.52516136106541430, 10.60460290274525023, 12.80182748008146961,
  15.10441257307551530, 17.50230784587388584, 19.98721449566188615,
  22.55216385312342289, 25.19122118273868150, 27.89927138384089157,
  30.67186010608067280, 33.50507345013688888, 36.39544520803305358,
  39.33988418719949404, 42.33561646075348503, 45.38013889847690803,
  48.47118135183522388, 51.60667556776437357, 54.78472939811231919,
  58.00360522298051994, 61.26170176100200198, 64.55753862700633106,
  67.88974313718153498, 71.25703896716800901];
phantasus.FisherExact.logFactorial = function (k) {
  if (k >= 30) { // stirlings approximation
    var C0 = 9.18938533204672742e-01;
    var C1 = 8.33333333333333333e-02;
    var C3 = -2.77777777777777778e-03;
    var C5 = 7.93650793650793651e-04;
    var C7 = -5.95238095238095238e-04;
    var r = 1.0 / k;
    var rr = r * r;
    return (k + 0.5) * Math.log(k) - k + C0 + r
      * (C1 + rr * (C3 + rr * (C5 + rr * C7)));
    // log k! = (k + 1/2)log(k) - k + (1/2)log(2Pi) + stirlingCorrection(k)
  }
  return phantasus.FisherExact.logFactorials[k];
};

phantasus.FisherExact.fisherTest = function (a, b, c, d) {
  // match R 2-sided fisher.test
  var p = phantasus.FisherExact.phyper(a, b, c, d);
  var sum = p;
  for (var _a = 0, n = a + b + c + d; _a <= n; _a++) {
    var _b = a + b - _a;
    var _c = a + c - _a;
    var _d = b + d - _b;
    if (_a !== a && _b >= 0 && _c >= 0 && _d >= 0) {
      var _p = phantasus.FisherExact.phyper(_a, _b, _c, _d);
      if (_p <= p) {
        sum += _p;
      }
    }
  }
  return Math.min(1, sum);
  // var lt = jStat.hypgeom.cdf(a, a + b + c + d, a + b, a + c);
  // var gt = jStat.hypgeom.cdf(b, a + b + c + d, a + b, b + d);
  // return Math.min(1, 2 * Math.min(lt, gt));
};
phantasus.FisherExact.toString = function () {
  return 'Fisher Exact Test';
};

phantasus.FoldChange = function (list1, list2) {
  var m1 = phantasus.Mean(list1);
  var m2 = phantasus.Mean(list2);
  return (m1 / m2);
};
phantasus.FoldChange.toString = function () {
  return 'Fold Change';
};
phantasus.TTest = function (list1, list2) {
  var m1 = phantasus.Mean(list1);
  var m2 = phantasus.Mean(list2);
  var s1 = Math.sqrt(phantasus.Variance(list1, m1));
  var s2 = Math.sqrt(phantasus.Variance(list2, m2));
  var n1 = phantasus.CountNonNaN(list1);
  var n2 = phantasus.CountNonNaN(list2);
  return ((m1 - m2) / Math.sqrt((s1 * s1 / n1) + (s2 * s2 / n2)));
};
phantasus.TTest.toString = function () {
  return 'T-Test';
};
phantasus.Spearman = function (list1, list2) {
  var flist1 = [];
  var flist2 = [];
  for (var i = 0, n = list1.size(); i < n; i++) {
    var v1 = list1.getValue(i);
    var v2 = list2.getValue(i);
    if (isNaN(v1) || isNaN(v2)) {
      continue;
    }
    flist1.push(v1);
    flist2.push(v2);
  }
  var rank1 = phantasus.Ranking(flist1);
  var rank2 = phantasus.Ranking(flist2);
  return phantasus.Pearson(new phantasus.Vector('', rank1.length)
    .setArray(rank1), new phantasus.Vector('', rank2.length)
    .setArray(rank2));
};
phantasus.Spearman.toString = function () {
  return 'Spearman rank correlation';
};
phantasus.WeightedMean = function (weights, values) {
  var numerator = 0;
  var denom = 0;
  for (var i = 0, size = values.size(); i < size; i++) {
    var value = values.getValue(i);
    if (!isNaN(value)) {
      var weight = Math.abs(weights.getValue(i));
      if (!isNaN(weight)) {
        numerator += (weight * value);
        denom += weight;
      }
    }
  }
  return denom === 0 ? NaN : numerator / denom;
};
phantasus.WeightedMean.toString = function () {
  return 'Weighted average';
};

phantasus.createOneMinusMatrixValues = function (dataset) {
  var f = function (listOne, listTwo) {
    return 1 - dataset.getValue(listOne.getIndex(), listTwo.getIndex());
  };
  f.toString = function () {
    return 'One minus matrix values (for a precomputed similarity matrix)';
  };
  return f;
};

phantasus.Pearson = function (listOne, listTwo) {
  var sumx = 0;
  var sumxx = 0;
  var sumy = 0;
  var sumyy = 0;
  var sumxy = 0;
  var N = 0;
  for (var i = 0, size = listOne.size(); i < size; i++) {
    var x = listOne.getValue(i);
    var y = listTwo.getValue(i);
    if (isNaN(x) || isNaN(y)) {
      continue;
    }
    sumx += x;
    sumxx += x * x;
    sumy += y;
    sumyy += y * y;
    sumxy += x * y;
    N++;
  }
  var numr = sumxy - (sumx * sumy / N);
  var denr = Math.sqrt((sumxx - (sumx * sumx / N))
    * (sumyy - (sumy * sumy / N)));
  return denr == 0 ? 1 : numr / denr;
};
phantasus.Pearson.toString = function () {
  return 'Pearson correlation';
};

phantasus.Jaccard = function (listOne, listTwo) {

  var orCount = 0;
  var andCount = 0;
  for (var i = 0, size = listOne.size(); i < size; i++) {
    var xval = listOne.getValue(i);
    var yval = listTwo.getValue(i);
    if (isNaN(xval) || isNaN(yval)) {
      continue;
    }
    var x = xval > 0;
    var y = yval > 0;
    if (x && y) {
      andCount++;
    } else if (x || y) {
      orCount++;
    }
  }
  if (orCount === 0) {
    return 1;
  }
  return 1 - (andCount / orCount);
};

phantasus.Jaccard.toString = function () {
  return 'Jaccard distance';
};

phantasus.Cosine = function (listOne, listTwo) {
  var sumX2 = 0;
  var sumY2 = 0;
  var sumXY = 0;
  for (var i = 0, size = listOne.size(); i < size; i++) {
    var x = listOne.getValue(i);
    var y = listTwo.getValue(i);
    if (isNaN(x) || isNaN(y)) {
      continue;
    }
    sumX2 += x * x;
    sumY2 += y * y;
    sumXY += x * y;
  }
  return (sumXY / Math.sqrt(sumX2 * sumY2));
};

phantasus.Cosine.toString = function () {
  return 'Cosine similarity';
};

phantasus.Euclidean = function (x, y) {
  var dist = 0;
  for (var i = 0, size = x.size(); i < size; ++i) {
    var x_i = x.getValue(i);
    var y_i = y.getValue(i);
    if (isNaN(x_i) || isNaN(y_i)) {
      continue;
    }
    dist += (x_i - y_i) * (x_i - y_i);
  }
  return Math.sqrt(dist);
};
phantasus.Euclidean.toString = function () {
  return 'Euclidean distance';
};
phantasus.OneMinusFunction = function (f) {
  var dist = function (x, y) {
    return 1 - f(x, y);
  };
  dist.toString = function () {
    var s = f.toString();
    return 'One minus ' + s[0].toLowerCase() + s.substring(1);
  };
  return dist;
};

phantasus.LinearRegression = function (xVector, yVector) {
  var sumX = 0;
  var sumY = 0;
  var sumXX = 0;
  var sumXY = 0;
  var count = 0;
  for (var i = 0, size = xVector.size(); i < size; i++) {
    var x = xVector.getValue(i);
    var y = yVector.getValue(i);
    if (!isNaN(x) && !isNaN(y)) {
      sumX += x;
      sumY += y;
      sumXX += x * x;
      sumXY += x * y;
      count++;
    }
  }

  var m = ((count * sumXY) - (sumX * sumY)) /
    ((count * sumXX) - (sumX * sumX));
  var b = (sumY / count) - ((m * sumX) / count);
  return {
    m: m,
    b: b
  };
};

phantasus.KendallsCorrelation = function (x, y) {

  /**
   * Returns the sum of the number from 1 .. n according to Gauss' summation formula:
   * \[ \sum\limits_{k=1}^n k = \frac{n(n + 1)}{2} \]
   *
   * @param n the summation end
   * @return the sum of the number from 1 to n
   */
  function sum(n) {
    return n * (n + 1) / 2;
  }

  var xArray = [];
  var yArray = [];
  for (var i = 0, size = x.size(); i < size; ++i) {
    var x_i = x.getValue(i);
    var y_i = y.getValue(i);
    if (isNaN(x_i) || isNaN(y_i)) {
      continue;
    }
    xArray.push(x_i);
    yArray.push(y_i);
  }
  var n = xArray.length;
  var numPairs = sum(n - 1);
  var pairs = [];
  for (var i = 0; i < n; i++) {
    pairs[i] = [xArray[i], yArray[i]];
  }
  pairs.sort(function (pair1, pair2) {
    var a = pair1[0];
    var b = pair2[0];
    var compareFirst = (a === b ? 0 : (a < b ? -1 : 1));
    if (compareFirst !== 0) {
      return compareFirst;
    }
    a = pair1[1];
    b = pair2[1];
    return (a === b ? 0 : (a < b ? -1 : 1));
  });

  var tiedXPairs = 0;
  var tiedXYPairs = 0;
  var consecutiveXTies = 1;
  var consecutiveXYTies = 1;
  var prev = pairs[0];
  for (var i = 1; i < n; i++) {
    var curr = pairs[i];
    if (curr[0] === prev[0]) {
      consecutiveXTies++;
      if (curr[1] === prev[1]) {
        consecutiveXYTies++;
      } else {
        tiedXYPairs += sum(consecutiveXYTies - 1);
        consecutiveXYTies = 1;
      }
    } else {
      tiedXPairs += sum(consecutiveXTies - 1);
      consecutiveXTies = 1;
      tiedXYPairs += sum(consecutiveXYTies - 1);
      consecutiveXYTies = 1;
    }
    prev = curr;
  }
  tiedXPairs += sum(consecutiveXTies - 1);
  tiedXYPairs += sum(consecutiveXYTies - 1);
  var swaps = 0;
  var pairsDestination = [];
  for (var segmentSize = 1; segmentSize < n; segmentSize <<= 1) {
    for (var offset = 0; offset < n; offset += 2 * segmentSize) {
      var i = offset;
      var iEnd = Math.min(i + segmentSize, n);
      var j = iEnd;
      var jEnd = Math.min(j + segmentSize, n);
      var copyLocation = offset;
      while (i < iEnd || j < jEnd) {
        if (i < iEnd) {
          if (j < jEnd) {
            var c = (pairs[i][1] === pairs[j][1] ? 0 : (pairs[i][1] < pairs[j][1] ? -1 : 1));
            if (c <= 0) {
              pairsDestination[copyLocation] = pairs[i];
              i++;
            } else {
              pairsDestination[copyLocation] = pairs[j];
              j++;
              swaps += iEnd - i;
            }
          } else {
            pairsDestination[copyLocation] = pairs[i];
            i++;
          }
        } else {
          pairsDestination[copyLocation] = pairs[j];
          j++;
        }
        copyLocation++;
      }
    }
    var pairsTemp = pairs;
    pairs = pairsDestination;
    pairsDestination = pairsTemp;
  }

  var tiedYPairs = 0;
  var consecutiveYTies = 1;
  prev = pairs[0];
  for (var i = 1; i < n; i++) {
    var curr = pairs[i];
    if (curr[1] === prev[1]) {
      consecutiveYTies++;
    } else {
      tiedYPairs += sum(consecutiveYTies - 1);
      consecutiveYTies = 1;
    }
    prev = curr;
  }
  tiedYPairs += sum(consecutiveYTies - 1);

  var concordantMinusDiscordant = numPairs - tiedXPairs - tiedYPairs + tiedXYPairs - 2 * swaps;
  var nonTiedPairsMultiplied = (numPairs - tiedXPairs) * (numPairs - tiedYPairs);
  return concordantMinusDiscordant / Math.sqrt(nonTiedPairsMultiplied);
};
phantasus.KendallsCorrelation.toString = function () {
  return 'Kendall\'s correlation';
};
