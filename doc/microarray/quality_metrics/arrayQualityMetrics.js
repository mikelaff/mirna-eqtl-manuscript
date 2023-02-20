// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = [ false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false ];
var arrayMetadata    = [ [ "1", "MA1.CEL", "E6", "2022-06-15", "D54-Control-3", "RNA30", "Control", "3", "Week 2", "2022-06-21", "1899-12-31 12:31:00", "1.97", "0.81", "0.038", "9.9", "2.3", "42.7" ], [ "2", "MA2.CEL", "E24", "2022-06-15", "D54-Control-6", "RNA48", "Control", "6", "Week 2", "2022-06-21", "1899-12-31 12:41:00", "2.05", "1.71", "0.018", "9.9", "2.4", "41.6" ], [ "3", "MA3.CEL", "E13", "2022-06-08", "D54-Control-2", "RNA13", "Control", "2", "Week 1", "2022-06-17", "1899-12-31 16:26:00", "2.25", "1.86", "0.033", "9.6", "2.2", "23.2" ], [ "4", "MA4.CEL", "E9", "2022-06-08", "D54-Control-1", "RNA9", "Control", "1", "Week 1", "2022-06-17", "1899-12-31 16:24:00", "2.15", "0.72", "0.01", "9.5", "2.2", "21.5" ], [ "5", "MA5.CEL", "E11", "2022-06-15", "D54-4707-6", "RNA35", "4707", "6", "Week 2", "2022-06-21", "1899-12-31 12:34:00", "2", "1.89", "0.093", "9.8", "2.3", "46.2" ], [ "6", "MA6.CEL", "E11", "2022-06-08", "D54-4707-4", "RNA11", "4707", "4", "Week 1", "2022-06-17", "1899-12-31 16:25:00", "2.05", "1.4", "0.07", "9.7", "2.3", "31.1" ], [ "7", "MA7.CEL", "E12", "2022-06-15", "D54-Control-4", "RNA36", "Control", "4", "Week 2", "2022-06-21", "1899-12-31 12:35:00", "2.05", "1.08", "0.025", "9.7", "2.2", "38.5" ], [ "8", "MA8.CEL", "E5", "2022-06-08", "D54-4707-2", "RNA5", "4707", "2", "Week 1", "2022-06-17", "1899-12-31 16:21:00", "2.09", "1.28", "0.026", "9.8", "2.3", "29.9" ], [ "9", "MA9.CEL", "E22", "2022-06-15", "D54-4707-2", "RNA46", "4707", "2", "Week 2", "2022-06-21", "1899-12-31 12:40:00", "2.02", "1.54", "0.012", "9.7", "2.3", "48.7" ], [ "10", "MA10.CEL", "E14", "2022-06-08", "D54-4707-5", "RNA14", "4707", "5", "Week 1", "2022-06-17", "1899-12-31 16:27:00", "2.28", "1.64", "0.022", "9.5", "2.2", "28.1" ], [ "11", "MA11.CEL", "E23", "2022-06-08", "D54-Control-6", "RNA23", "Control", "6", "Week 1", "2022-06-17", "1899-12-31 16:32:00", "2.2", "0.47", "0.019", "9.4", "2.1", "18.3" ], [ "12", "MA12.CEL", "E18", "2022-06-15", "D54-4707-4", "RNA42", "4707", "4", "Week 2", "2022-06-21", "1899-12-31 12:38:00", "2.07", "1.71", "0.015", "9.9", "2.6", "49.2" ], [ "13", "MA13.CEL", "E16", "2022-06-15", "D54-4707-5", "RNA40", "4707", "5", "Week 2", "2022-06-21", "1899-12-31 12:37:00", "2", "1.66", "0.011", "9.8", "2.4", "49" ], [ "14", "MA14.CEL", "E3", "2022-06-15", "D54-Control-5", "RNA27", "Control", "5", "Week 2", "2022-06-21", "1899-12-31 12:29:00", "2", "1.08", "0.054", "9.9", "2.5", "41.6" ], [ "15", "MA15.CEL", "E2", "2022-06-08", "D54-4707-1", "RNA2", "4707", "1", "Week 1", "2022-06-17", "1899-12-31 16:19:00", "2.18", "0.88", "0.041", "9.8", "2.2", "29.5" ], [ "16", "MA16.CEL", "E18", "2022-06-08", "D54-Control-4", "RNA18", "Control", "4", "Week 1", "2022-06-17", "1899-12-31 16:29:00", "2.16", "1.44", "0.016", "9.7", "2.2", "24" ], [ "17", "MA17.CEL", "E10", "2022-06-08", "D54-4707-3", "RNA10", "4707", "3", "Week 1", "2022-06-17", "1899-12-31 16:24:00", "2.09", "1.31", "0.034", "9.9", "2.3", "31.6" ], [ "18", "MA18.CEL", "E19", "2022-06-08", "D54-4707-6", "RNA19", "4707", "6", "Week 1", "2022-06-17", "1899-12-31 16:30:00", "2.2", "1.87", "0.031", "9.8", "2.2", "32.7" ], [ "19", "MA19.CEL", "E4", "2022-06-15", "D54-4707-1", "RNA28", "4707", "1", "Week 2", "2022-06-21", "1899-12-31 12:30:00", "2.08", "1.4", "0.015", "9.9", "2.5", "49.7" ], [ "20", "MA20.CEL", "E8", "2022-06-15", "D54-4707-3", "RNA32", "4707", "3", "Week 2", "2022-06-21", "1899-12-31 12:32:00", "2.06", "1.81", "0.009", "9.9", "2.5", "49.8" ], [ "21", "MA21.CEL", "E16", "2022-06-08", "D54-Control-3", "RNA16", "Control", "3", "Week 1", "2022-06-17", "1899-12-31 16:28:00", "2.11", "1.01", "0.038", "9.6", "2.4", "20.6" ], [ "22", "MA22.CEL", "E13", "2022-06-15", "D54-Control-1", "RNA37", "Control", "1", "Week 2", "2022-06-21", "1899-12-31 12:35:00", "2.07", "1.95", "0.024", "10", "2.5", "43.1" ], [ "23", "MA23.CEL", "E21", "2022-06-15", "D54-Control-2", "RNA45", "Control", "2", "Week 2", "2022-06-21", "1899-12-31 12:40:00", "2.02", "1.25", "0.015", "10", "2.3", "39.2" ], [ "24", "MA24.CEL", "E21", "2022-06-08", "D54-Control-5", "RNA21", "Control", "5", "Week 1", "2022-06-17", "1899-12-31 16:31:00", "2.27", "1.87", "0.058", "9.8", "2.1", "23.7" ] ];
var svgObjectNames   = [ "pca", "dens" ];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
    for(i=0; i<ssrules.length; i++) {
        if (ssrules[i].selectorText == (".aqm" + reportObjId)) {
		ssrules[i].style.cssText = cssText[0+status];
		break;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
