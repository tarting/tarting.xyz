// Modified from the example at https://medium.com/@mariusc23/hide-header-on-scroll-down-show-on-scroll-up-67bbaae9a78c
// Hide Header on on scroll down
var didScroll;
var didResize = true;
var lastScrollTop = 0;
var delta = 5;
var headerHeight;


$(window).scroll(function(event){
    didScroll = true;
});


setInterval(function() {
    if (didScroll) {
        hasScrolled();
        didScroll = false;
    }
}, 250);


function hasScrolled() {
    var st = $(this).scrollTop();
    if (Math.abs(lastScrollTop - st) <= delta) {
        lastScrollTop = st;
        return;
    } else if ((st-lastScrollTop) < 0) {
        $("header").css("top", 0);
    } else if (st <= $('header').outerHeight()/3.0) {
        $("header").css("top", 0);
    } else {
        $("header").css("top", -$("header").outerHeight());
    }

    lastScrollTop = st;
    return;
}

// https://stackoverflow.com/a/14254924
function jq( myid ) {
    return myid.replace( /(:|\.|\[|\]|,|=|@)/g, "\\$1" );
}

$(document).on('click', 'a[href^="#"]', function(e) {
//$("footnote-ref").click(function() {
    e.preventDefault();
    var id = jq($(this).attr("href"));
    var link = $(this);
    var target = $(id);
    console.log(link.text);

    // https://stackoverflow.com/a/42769683
    var fontSize = parseFloat(getComputedStyle(document.documentElement).fontSize);

    var off_top = target.offset().top - $('header').outerHeight() - fontSize;
    $('body, html').animate({scrollTop: off_top});
});
