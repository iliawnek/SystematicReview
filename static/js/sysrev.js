$(function () {
    $('[data-toggle="tooltip"]').tooltip();
});

function toggleChevron(e) {
    $(".collapse-button").toggleClass('glyphicon-chevron-down glyphicon-chevron-up');
}
$('.collapse').on('hidden.bs.collapse', toggleChevron);
$('.collapse').on('shown.bs.collapse', toggleChevron);
