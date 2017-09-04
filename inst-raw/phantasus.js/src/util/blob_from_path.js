phantasus.BlobFromPath = function () {
};
phantasus.BlobFromPath.getFileBlob = function (url, cb) {
  var xhr = new XMLHttpRequest();
  xhr.open("GET", url);
  xhr.responseType = "blob";
  xhr.addEventListener('load', function () {
    cb(xhr.response);
  });
  xhr.send();
};

phantasus.BlobFromPath.blobToFile = function (blob) {
  blob.lastModifiedDate = new Date();
  return blob;
};

phantasus.BlobFromPath.getFileObject = function (filePathOrUrl, cb) {
  phantasus.BlobFromPath.getFileBlob(filePathOrUrl, function (blob) {
    cb(phantasus.BlobFromPath.blobToFile(blob));
  });
};