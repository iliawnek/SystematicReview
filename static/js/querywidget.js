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
                operators: ['contains', 'not_contains']
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

    function getRulesFromQuery(query) {
        var rules = {};
        return rules;
    }

    function getQueryFromRules(outer) {
        if (!outer)
            return undefined;

        var condition = outer.condition;
        if (!condition) {
            console.error("No condition for");
            console.error(outer);
            return undefined;
        }

        outer = outer.rules;

        var parts = [];
        for (var i = 0; i < outer.length; i++) {
            var rule = outer[i];

            console.log("Adding " + rule);

            if (rule.condition) {
                parts.push('(' + getQueryFromRules(rule) + ')');
                continue;
            }

            var value = rule.value;
            var field = rule.field;

            switch (rule.operator) {
                case 'contains':
                    break;
                case 'not_contains':
                    value = "NOT " + value;
                    break;
                default:
                    console.warn("Unsupported operator: " + rule.operator);
            }

            console.log(value);
            parts.push(value + (field == 'All Fields' ? '' : '[' + field + ']'));
        }

        console.log(parts);
        return parts.join(' ' + condition + ' ');
    }

    function getQuery(it) {
        var rules = it.queryBuilder('getRules');
        console.log(rules);
        return getQueryFromRules(rules);
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

        var textArea = it;
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

        var timeout = false;
        var update = function () {
            if (timeout) {
                clearTimeout(timeout);
                timeout = false;
            }

            timeout = setTimeout(function () {
                var focused = $(':focus');

                if (focused[0] && focused[0].tagName == 'SELECT')
                    return;

                // workaround for getRules not working. not sure why
                textArea.focus();

                var query = getQuery(it);
                console.log(query);
                if (query)
                    textArea.val(query);

                focused.focus();
            }, 1000);
        };
        it.keydown(update);
        it.click(update);
    });
});
