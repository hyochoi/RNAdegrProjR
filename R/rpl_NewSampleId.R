#' Replace colnames of a matrix based on NewSampleId in Pilot2
#'
#' @param mat a matrix with entity_submitter_id columns
#' @param Pilot2 a reference matrix having NewSampleId
#' @export

rpl_NewSampleId = function(mat, Pilot2) {
  tmat <- t(mat)
  OldSampleId <- rownames(tmat)
  mat0 <- cbind(OldSampleId, tmat)

  # Keep the selected patients in Pilot data
  OldNewId <- Pilot2[, c("entity_submitter_id", "NewSampleId")]
  mat1 <- merge(OldNewId, mat0, by.x="entity_submitter_id", by.y="OldSampleId")

  # Change to a NewSampleId, sort by a NewSampleId
  mat2 <- mat1[order(mat1$NewSampleId), ]
  mat3 <- mat2[, 3:dim(mat2)[2]]
  rownames(mat3) <- c(mat2$NewSampleId)
  mat4 <- t(mat3)
  mat5 <- matrix(as.numeric(mat4), ncol=ncol(mat4))
  rownames(mat5) <- row.names(mat4)
  colnames(mat5) <- colnames(mat4)

  return(mat5)
}
