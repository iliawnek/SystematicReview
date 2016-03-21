$(function () {
    $('[data-toggle="tooltip"]').tooltip();
});

$(function () {
    $(".collapse-button").click(function () {
        $(".collapse-button").toggleClass("glyphicon-menu-down glyphicon-menu-up")
    })
});