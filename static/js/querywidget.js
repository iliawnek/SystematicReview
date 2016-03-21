// see sysrev.widgets.QueryWidget

$(function () {
    if (window.queryWidgetIncluded)
        return;

    // To dump fields from PubMed site:
    // var options = document.querySelectorAll("#ff_0 > option");
    // var res = [];
    // for (var i = 0; i < options.length; i++) {
    //     res.push(options[i].getAttribute("value"));
    // }
    // console.log("'" + res.join("', '") + "'");

    var filters = [
        'Affiliation', 'All Fields', 'Author', 'Author - Corporate', 'Author - First', 'Author - Full',
        'Author - Identifier', 'Author - Last', 'Book', 'Date - Completion', 'Date - Create', 'Date - Entrez',
        'Date - MeSH', 'Date - Modification', 'Date - Publication', 'EC/RN Number', 'Editor', 'Filter',
        'Grant Number', 'ISBN', 'Investigator', 'Investigator - Full', 'Issue', 'Journal', 'Language',
        'Location ID', 'MeSH Major Topic', 'MeSH Subheading', 'MeSH Terms', 'Other Term', 'Pagination',
        'Pharmacological Action', 'Publication Type', 'Publisher', 'Secondary Source ID', 'Subject - Personal Name',
        'Supplementary Concept', 'Text Word', 'Title', 'Title/Abstract', 'Transliterated Title', 'Volume'
    ];

    var convertedFilters;

    function buildFilters() {
        if (convertedFilters)
            return convertedFilters;

        convertedFilters = [];

        $.each(filters, function (i, it) {
            convertedFilters.push({
                id: it,
                label: it,
                type: 'string',
                operators: ['contains']
            });
        });

        return convertedFilters;
    }

    function displayBuilder(it) {
        it.queryBuilder({
            // TODO bootstrap tooltip error plugin? requires Bootstrap Tooltip
            // plugins: ['bt-tooltip-errors'],

            filters: buildFilters()

            // TODO: validation rules?
            // rules: rules_basic
        });
    }

    window.queryWidgetIncluded = true;

    var uid = 0;
    $("textarea.queryWidget").each(function (i, it) {
        it = $(it);
        var id = uid++;
        var advanced = false;
        console.log("Found queryWidget textarea");
        console.log(it);

        var qwId = "queryWidget_" + id;

        it.before($("<div id='" + qwId + "' class='queryWidget'></div>"));

        it = $('#' + qwId);

        function toggleAdvanced() {
            advanced = !advanced;
            if (advanced) {
                displayBuilder(it);
            } else {
                it.clear();
            }
        }

        toggleAdvanced();

        $('#target_' + qwId).droppable({
            drop: function (e, ui) {
                console.log("Dropped %o %o", e, ui);
                $(ui.draggable).clone().appendTo($(this));
            }
        });
    });
});
