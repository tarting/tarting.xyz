// Modified from the example at https://medium.com/@mariusc23/hide-header-on-scroll-down-show-on-scroll-up-67bbaae9a78c
// Hide Header on on scroll down
var didScroll;
var lastScrollTop = 0;
var delta = 5;
var navbarHeight = $('header').outerHeight();

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
    console.log(st, lastScrollTop, st-lastScrollTop);
    if (Math.abs(lastScrollTop - st) <= delta) {
        lastScrollTop = st;
        return;
    } else if ((st-lastScrollTop) < 0) {
        $('header').removeClass('nav-up').addClass('nav-down');
    } else if (st <= $('header').outerHeight()/3.0) {
        $('header').removeClass('nav-up').addClass('nav-down');
    } else {
        $('header').removeClass('nav-down').addClass('nav-up');
    }

    lastScrollTop = st;
    return;
}
