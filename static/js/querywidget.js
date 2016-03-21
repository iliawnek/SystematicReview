// see sysrev.widgets.QueryWidget

$(function() {
    if (window.queryWidgetIncluded)
        return;

    var types = ['AND', 'OR'];

    function addWidgets(it) {
        $.each(types, function(i, type) {
            // TODO: Add styles so these make more sense
            // TODO: Add destination spots on side for text/other widgets?
            it.append($("<div class='queryWidgetType'><span class='queryWidgetTarget queryWidgetTargetLeft'></span>" + type + "<span class='queryWidgetTarget queryWidgetTargetRight'></div>"))
        });
        it.find('.queryWidgetType').draggable({
            revert:"valid",
            helper:"clone"
        });
    }

    window.queryWidgetIncluded = true;

    var uid = 0;
    $("textarea.queryWidget").each(function(i, it) {
        it = $(it);
        var id = uid++;
        var advanced = false;
        console.log("Found queryWidget textarea");
        console.log(it);

        var qwId = "queryWidget_" + id;

        it.before($("<div id='" + qwId + "' class='queryWidget'></div>"));
        //TODO: remove inline test style
        it.before($("<div id='target_" + qwId + "' class='queryWidgetTarget' style='background-color: red; min-width: 100px; min-height: 10px;'></div>"));

        it = $('#' + qwId);

        addWidgets(it);

        $('#target_' + qwId).droppable({
            drop: function (e, ui) {
                console.log("Dropped %o %o", e, ui);
                $(ui.draggable).clone().appendTo($(this));
            }
        });
    });
});
