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