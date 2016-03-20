$(function () {
  $('[data-toggle="tooltip"]').tooltip();
});

$(function () {
  $("#reading").click(function(){
    $(".panel-body").toggleClass("panel-reading");
    $("#reading").toggleClass("active");
  });
  $('[data-toggle="tooltip"]').tooltip();
});