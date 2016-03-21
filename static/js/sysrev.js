$(function () {
    $('[data-toggle="tooltip"]').tooltip();
});

$(function () {
    $(".collapse-button").click(function () {
        $(".collapse-button").toggleClass("glyphicon-chevron-down glyphicon-chevron-up")
    })
});