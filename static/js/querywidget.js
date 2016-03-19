// see sysrev.widgets.QueryWidget

$(function() {
    if (window.queryWidgetIncluded)
        return;

    window.queryWidgetIncluded = true;

    $("textarea.queryWidget").each(function(i, it) {
        console.log("Found queryWidget textarea");
        console.log(it);
    });
});