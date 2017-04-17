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

function checkUploadFiles(target, cssClass, fileExtension, maxSize) {
	if (target[0].files[0].size > maxSize) {
		showMessage('The ' + cssClass + ' file is too large.', cssClass);
	} else if (!endsWith(target.val().toLowerCase(), fileExtension)) {
		console.log(target.val().toLowerCase());
		showMessage('The ' + cssClass + ' file must be a ' + fileExtension + ' file.', cssClass);
	} else {
		clearMessages(cssClass);
	}
}

$('[name=alignment]').bind('change', function() {
	target = $('[name=alignment]');
	maxSize = 16 * 1024 * 1024; // 20MB
	checkUploadFiles(target, 'alignment', '.sam', maxSize)
});

$('[name=reference]').bind('change', function() {
	target = $('[name=reference]');
	maxSize = 3 * 1024 * 1024; // 3MB
	checkUploadFiles(target, 'reference', '.csv', maxSize)
});

$(function() {
    $(".accordion").accordion({
        animate: false,
        collapsible: true,
        heightStyle: "content",
    });
});
