phantasus.VectorKeys = {};
/** [string] of field names in array */
phantasus.VectorKeys.FIELDS = 'phantasus.fields';
phantasus.VectorKeys.VALUE_TO_INDICES = 'phantasus.valueToIndices';
/** [int] of visible field indices in phantasus.VectorKeys.FIELDS */
phantasus.VectorKeys.VISIBLE_FIELDS = 'phantasus.visibleFields';
phantasus.VectorKeys.DATA_TYPE = 'phantasus.dataType';
/** Function to map an array to a single value for sorting */
phantasus.VectorKeys.ARRAY_SUMMARY_FUNCTION = 'phantasus.arraySummaryFunct';
/** Key for object (e.g. box plot) that summarizes data values */
phantasus.VectorKeys.HEADER_SUMMARY = 'phantasus.headerSummary';
/** Key indicating to show header summary */
phantasus.VectorKeys.SHOW_HEADER_SUMMARY = 'phantasus.showHeaderSummary';

phantasus.VectorKeys.TITLE = 'phantasus.title';
/** Function to compute vector value */
phantasus.VectorKeys.FUNCTION = 'phantasus.funct';

/** Indicates that vector values are dynamically computed based on selection */
phantasus.VectorKeys.SELECTION = 'phantasus.selection';

/** Whether to recompute a function when creating a new heat map */
phantasus.VectorKeys.RECOMPUTE_FUNCTION = 'phantasus.recompute.funct';

/** Boolean, whether to recompute a function when heat map selection changes */
phantasus.VectorKeys.RECOMPUTE_FUNCTION_SELECTION = 'phantasus.recompute.funct.selection';

phantasus.VectorKeys.COPY_IGNORE = new phantasus.Set();
phantasus.VectorKeys.COPY_IGNORE.add(phantasus.VectorKeys.HEADER_SUMMARY);
phantasus.VectorKeys.COPY_IGNORE.add(phantasus.VectorKeys.DATA_TYPE);
phantasus.VectorKeys.COPY_IGNORE.add(phantasus.VectorKeys.VALUE_TO_INDICES);

phantasus.VectorKeys.JSON_WHITELIST = new phantasus.Set();
phantasus.VectorKeys.JSON_WHITELIST.add(phantasus.VectorKeys.FIELDS);
phantasus.VectorKeys.JSON_WHITELIST.add(phantasus.VectorKeys.DATA_TYPE);
phantasus.VectorKeys.JSON_WHITELIST.add(phantasus.VectorKeys.TITLE);
