MAX_ALIGNMENT_FILE_SIZE = 80 * 1024 * 1024; // 80MB
MAX_REFERENCE_FILE_SIZE = 5 * 1024 * 1024; // 5MB

function endsWith(str, suffix) {
    return str.indexOf(suffix, str.length - suffix.length) !== -1;
}

function showMessage(message, cssClass) {
    clearMessages(cssClass);
    $(".flashes").append("<ul class='" + cssClass + "'><li>" + message + "</li></ul>");
}

function clearMessages(cssClass) {
    $("." + cssClass).remove();
}

function checkUploadFiles(target, description, fileExtension, maxSize, cssClass) {
    clearMessages("server-error");
    if (target[0].files[0].size > maxSize) {
        showMessage('The ' + description + ' file is too large.', cssClass);
    } else if (!endsWith(target.val().toLowerCase(), fileExtension)) {
        showMessage('The ' + description + ' file must be a ' + fileExtension + ' file.', cssClass);
    } else {
        clearMessages(cssClass);
    }
}

$(function() {
    $(".accordion").accordion({
        animate: false,
        collapsible: true,
        heightStyle: "content",
    });
});

$(document).ready(function() {
    $('[name=alignment]').bind('change', function() {
        target = $('[name=alignment]');
        checkUploadFiles(target, 'alignment', '.sam', MAX_ALIGNMENT_FILE_SIZE, 'alignment-file-error');
    });

    $('[name=reference]').bind('change', function() {
        target = $('[name=reference]');
        checkUploadFiles(target, 'reference', '.csv', MAX_REFERENCE_FILE_SIZE, 'reference-file-error');
    });
});
