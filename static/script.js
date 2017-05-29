MAX_ALIGNMENT_FILE_SIZE = 40 * 1024 * 1024; // 40MB
MAX_REFERENCE_FILE_SIZE = 5 * 1024 * 1024; // 5MB

PERMISSION_TO_UPLOAD_REFERENCE = false;
PERMISSION_TO_UPLOAD_ALIGNMENT = false;

function endsWith(str, suffix) {
    return str.indexOf(suffix, str.length - suffix.length) !== -1;
}

function showMessage(message, cssClass) {
    $(".flashes").append("<ul class='" + cssClass + "'><li>" + message + "</li></ul>");
}

function clearMessages(cssClass) {
    $("." + cssClass).remove();
}

function checkUploadFiles(target, description, fileExtension, maxSize, cssClass, permission) {
    clearMessages("server-error");
	fileOfInterest = target[0].files[0]
    if (fileOfInterest.size > maxSize) {
		clearMessages(cssClass);
        showMessage('The ' + description + ' file is too large.', cssClass);
    } else if (!endsWith(target.val().toLowerCase(), fileExtension)) {
		clearMessages(cssClass);
        showMessage('The ' + description + ' file must be a ' + fileExtension + ' file.', cssClass);
    } else {
		clearMessages(cssClass);
    }
}

function getSignedRequest(file, description){
	if (!file) {
		showMessage("No " + description + " file.", "server-error");
	} else {
		$.get("/sign_s3/", {"type": file.type})
			.fail(function() {
				showMessage("Server error: failed to get signed request for " + description + ".", "server-error");
			})
			.done(function (data) {
				parsedData = JSON.parse(data)
				uploadFile(file, parsedData.data, parsedData.url, description);
			});
	}
}

function setHiddenFieldValues(field, value){
	$("[name=" + field + "-hidden]")[0].value = value;
}

function uploadFile(file, s3, url, description){
	console.log(s3);

	s3Form = new FormData();
	for(key in s3.fields){
		s3Form.append(key, s3.fields[key]);
	}
	s3Form.append('file', file);

	// $.post(s3.url, s3Form)
	// 	.fail(function (){
	// 		showMessage("Failed to upload " + description + " file.", "server-error");
	// 	})
	// 	.done(function (data){
	// 		setHiddenFieldValues(description, url)
	// 	});

	$.ajax({
        url: s3.url,
        type: 'POST',

        data: s3Form,

        cache: false,
        contentType: false,
        processData: false,
    })
		.fail(function (){
			showMessage("Failed to upload " + description + " file.", "server-error");
		})
		.done(function (data){
			setHiddenFieldValues(description, url)
		});
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
        checkUploadFiles(target, 'alignment', '.sam', MAX_ALIGNMENT_FILE_SIZE, 'alignment-file-error', PERMISSION_TO_UPLOAD_ALIGNMENT);
    });

    $('[name=reference]').bind('change', function() {
        target = $('[name=reference]');
        checkUploadFiles(target, 'reference', '.csv', MAX_REFERENCE_FILE_SIZE, 'reference-file-error', PERMISSION_TO_UPLOAD_REFERENCE);
    });

	$('[name=submit]').on('click', function(e) {
		clearMessages("server-error");
		e.preventDefault();

		$.when(getSignedRequest($('[name=alignment]')[0].files[0], 'alignment', e),
			   getSignedRequest($('[name=reference]')[0].files[0], 'reference', e))
			.done(function() {
				if ($("[name=reference-hidden]")[0].value != "" && $("[name=alignment-hidden]")[0].value != ""){
					$('form').submit();
				}
			});
	});
});
