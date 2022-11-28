$(document).ready(function () {

	$("body").tooltip({ selector: '[data-toggle=tooltip]' });
      
    $("[data-toggle=tooltip]").tooltip({html:true});
      
    $("label").each(function () {
      	$(this).html($(this).html().replaceAll("&lt;", "<").replaceAll("&gt;", ">"));
    });
      
	$('.moveTo').click(function (event) {
		var tab = $(this).attr('href');
		$('#helptabs').find('a[href="'+tab+'"]').click()
	});
      
	$('.moveToRef').click(function (event) {
   		var ref = $(this).attr('href');
       	$('#helptabs').find('a[href="'+"#refs"+'"]').click();
       	$(ref).css("background-color", "#E4F1EF").fadeIn("slow");
		window.setTimeout(function() {
       		/*$(ref).removeAttr('style');*/
       		$(ref).css("background-color", "").fadeIn("slow");
   		}, 1000);
   	});
   	
    $('.gallery_pics').click(function(e) {
    	$(this).toggleClass('fullscreen');
  	});
  	
  	document.getElementById("CountsRate-setCV2").parentNode.appendChild(document.getElementById("CountsRate-CV"));
  	document.getElementById("CountsRate-CV").style.cssText+="display: inline-block; width: auto;";
  	document.getElementById("CountsRate-setCV2").nextElementSibling.style.cssText+="top:9px;";
  	
  	document.getElementById("CountsStrain1-setCV2").parentNode.appendChild(document.getElementById("CountsStrain1-CV"));
  	document.getElementById("CountsStrain1-CV").style.cssText+="display: inline-block; width: auto;";
  	/*document.getElementById("CountsStrain1-setCV2").nextElementSibling.style.cssText+="top:9px;";*/
  	
  	document.getElementById("CountsStrain2-setCV2").parentNode.appendChild(document.getElementById("CountsStrain2-CV"));
  	document.getElementById("CountsStrain2-CV").style.cssText+="display: inline-block; width: auto;";
  	/*document.getElementById("CountsStrain2-setCV2").nextElementSibling.style.cssText+="top:9px;";*/
      
});

$(window).on("load", function(){
    $(".selectize-input input").attr('readonly','readonly');
});
      
$('a[data-toggle="tab"]').on( 'shown.bs.tab', function() {
    $(window).resize();
});

window.FontAwesomeConfig = {
    searchPseudoElements: true
};
  	
function ToggleClassFold() {
	document.querySelectorAll("div[data-value='Strain 1']").forEach(element=>{element.className = "btn btn-warning ui-sortable-handle"});
	document.querySelectorAll("div[data-value='Strain 2']").forEach(element=>{element.className = "btn btn-warning ui-sortable-handle"});
	document.querySelectorAll("div[data-value='Strain 3']").forEach(element=>{element.className = "btn btn-warning ui-sortable-handle"});
	document.querySelectorAll("div[data-value='Strain 4']").forEach(element=>{element.className = "btn btn-warning ui-sortable-handle"});
	document.querySelectorAll("div[data-value='Strain 5']").forEach(element=>{element.className = "btn btn-warning ui-sortable-handle"});
	document.querySelectorAll("div[data-value='Strain 6']").forEach(element=>{element.className = "btn btn-warning ui-sortable-handle"});
};

function FoldHideFeedback() {
	const eq = document.getElementById("FoldEquation");
	var el1 = document.getElementById("FoldEquation-icon");
	if (el1 != null) el1.parentElement.removeChild(el1);
	
	var el2 = document.getElementById("FoldEquation-text");
	if (el2 != null) el2.parentElement.removeChild(el2);
	
	eq.style.border = '1px solid #cccccc';
	
	document.getElementById("FoldEquation-label").style.color = '';
// 	eq.style.paddingRight = '';

};

function FoldShowFeedback(text) {
	const eq = document.getElementById("FoldEquation");
	FoldHideFeedback()
	
	document.getElementById("FoldEquation-label").style.color = 'rgb(185, 74, 72)';
	
	let el1 = document.createElement('span');
	el1.className = 'form-control-feedback';
	el1.setAttribute('id', 'FoldEquation-icon');
	el1.setAttribute('style', 'color: #b94a48; right: 15px; top: 25px;');

	let el2 = document.createElement('i');
	el2.className = 'fas fa-triangle-exclamation';
	el2.setAttribute('role', 'presentation');
	el2.setAttribute('aria-label', 'triangle-exclamation icon');

	el1.appendChild(el2);

	eq.parentElement.appendChild(el1);
	
	let el3 = document.createElement('div');
	el3.setAttribute('id', 'FoldEquation-text');
	
	let el4 = document.createElement('p');
	el4.setAttribute('style', 'color: #b94a48; margin-top: 0px;');
	
	let el5 = document.createTextNode(text);
	
	el4.appendChild(el5);
	el3.appendChild(el4);
	
	eq.parentElement.appendChild(el3);
	
	eq.style.border = '1px solid rgb(185, 74, 72)';
// 	eq.style.paddingRight = '47.5px';

};

function FoldRemoveElements() {	
	const eq = document.getElementById("FoldEquation");
  	while (eq.firstChild) {
    	eq.removeChild(eq.lastChild);
  	};
};

function SwitchTabToPower() {
	$('a[data-value="Power"]').click();
};
