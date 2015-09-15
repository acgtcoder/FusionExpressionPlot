context("systemRunUtils.R")

describe("is.safeFileName()", {

   alwaysGoodNames <- c("x", "a.txt", "3", "/me/and/you", "this_or-that")
   alwaysBadNames <- c("`ls`", "$(USER)", "$", "(", ")", "f.txt;ls", "|ls",
      "> x.txt", "<bad.txt", "&ls", "+opt", "-opt", "arg +opt", "arg -opt")
   defaultGoodNames <- c("blank ok", "any~ok" )
   defaultBadNames <- c("no,", "no:", "no\\", "no\\\\")
   okWeakName <- "OKweak%*@'=[]{}?!"
   existingFile <- tempfile("deleteMeAfterTesting");
   file.create(existingFile);

   describe("default behavior (white listing)", {
      it("accepts always good characters", {
         got <- is.safeFileName( alwaysGoodNames )
         expect_true( all(got) )
      })
      it("accepts default good characters", {
         got <- is.safeFileName( defaultGoodNames )
         expect_true( all(got) )
      })
      it("rejects always bad characters", {
         got <- is.safeFileName( alwaysBadNames )
         expect_false( all(got) )
      })
      it("rejects default bad characters", {
         got <- is.safeFileName( defaultBadNames )
         expect_false( all(got) )
      })
      it("rejects weakly safe names", {
         got <- is.safeFileName( okWeakName )
         expect_false(got)
      })
      it("allows both existing and non existing file names", {
         expect_true( file.exists (existingFile ))
         got <- is.safeFileName( existingFile )
         expect_true(got)
         expect_false( file.exists("noSuchFile.exists" ))
         got <- is.safeFileName( "noSuchFile.exists" )
         expect_true(got)
      })

   }) # END: default

   describe("weak= FALSE white-listing (explicit default)", {
      describe("selectable options - default", {
         it("Accepts always good characters", {
            got <- is.safeFileName( alwaysGoodNames, weak= FALSE )
            expect_true( all(got) )
         })
         it("accepts default good characters", {
            got <- is.safeFileName( defaultGoodNames, weak= FALSE )
            expect_true( all(got) )
         })
         it("rejects always bad characters", {
            got <- is.safeFileName( alwaysBadNames, weak= FALSE )
            expect_false( all(got) )
         })
         it("rejects default bad characters", {
            got <- is.safeFileName( defaultBadNames, weak= FALSE )
            expect_false( all(got) )
         })
         it("rejects weakly safe names", {
            got <- is.safeFileName( okWeakName, weak= FALSE )
            expect_false(got)
         })
      }) # END: default selectable
      describe("selectable options - explicit defaults", {
         it("Accepts always good characters", {
            got <- is.safeFileName( alwaysGoodNames, weak= FALSE,
               okSpace= TRUE, okTilde= TRUE, okBackslash=FALSE, okComma=FALSE, okColon=FALSE)
            expect_true( all(got) )
         })
         it("accepts default good characters", {
            got <- is.safeFileName( defaultGoodNames, weak= FALSE,
               okSpace= TRUE, okTilde= TRUE, okBackslash=FALSE, okComma=FALSE, okColon=FALSE)
            expect_true( all(got) )
         })
         it("rejects always bad characters", {
            got <- is.safeFileName( alwaysBadNames, weak= FALSE,
               okSpace= TRUE, okTilde= TRUE, okBackslash=FALSE, okComma=FALSE, okColon=FALSE)
            expect_false( all(got) )
         })
         it("rejects default bad characters", {
            got <- is.safeFileName( defaultBadNames, weak= FALSE,
               okSpace= TRUE, okTilde= TRUE, okBackslash=FALSE, okComma=FALSE, okColon=FALSE)
            expect_false( all(got) )
         })
         it("rejects weakly safe names", {
            got <- is.safeFileName( okWeakName, weak= FALSE,
               okSpace= TRUE, okTilde= TRUE, okBackslash=FALSE, okComma=FALSE, okColon=FALSE)
            expect_false(got)
         })
      }) # END: selectable explicit defaults
      describe("selectable options - reverse of defaults", {
         it("Accepts always good characters", {
            got <- is.safeFileName( alwaysGoodNames, weak= FALSE,
               okSpace= FALSE, okTilde= FALSE, okBackslash=TRUE, okComma=TRUE, okColon=TRUE)
            expect_true( all(got) )
         })
         it("accepts default good characters", {
            got <- is.safeFileName( defaultGoodNames, weak= FALSE,
               okSpace= FALSE, okTilde= FALSE, okBackslash=TRUE, okComma=TRUE, okColon=TRUE)
            expect_false( all(got) )
         })
         it("rejects always bad characters", {
            got <- is.safeFileName( alwaysBadNames, weak= FALSE,
               okSpace= FALSE, okTilde= FALSE, okBackslash=TRUE, okComma=TRUE, okColon=TRUE);
            expect_false( all(got) )
         })
         it("rejects default bad characters", {
            got <- is.safeFileName( defaultBadNames, weak= FALSE,
               okSpace= FALSE, okTilde= FALSE, okBackslash=TRUE, okComma=TRUE, okColon=TRUE)
            expect_true( all(got) )
         })
         it("rejects weakly safe names", {
            got <- is.safeFileName( okWeakName, weak= FALSE,
               okSpace= FALSE, okTilde= FALSE, okBackslash=TRUE, okComma=TRUE, okColon=TRUE)
            expect_false(got)
         })
      }) # END: selectable reversed from defaults
      describe("Exists option in context of weak=FALSE",{
         expect_true( file.exists (existingFile ))
         expect_false( file.exists("noSuchFile.exists" ))

         it("doesn't affect default exists option", {
            got <- is.safeFileName( existingFile, weak= FALSE )
            expect_true(got)
            got <- is.safeFileName( "noSuchFile.exists", weak= FALSE)
            expect_true(got)
         })
         it("doesn't affect exists = TRUE", {
            got <- is.safeFileName( existingFile, exists= TRUE, weak= FALSE )
            expect_true(got)
            got <- is.safeFileName( "noSuchFile.exists", exists= TRUE, weak= FALSE)
            expect_false(got)
         })
         it("doesn't affect exists = FALSE", {
            got <- is.safeFileName( existingFile, exists= FALSE, weak= FALSE )
            expect_false(got)
            got <- is.safeFileName( "noSuchFile.exists", exists= FALSE, weak= FALSE)
            expect_true(got)
         })
      }) # END exist=
   }) # END weak=FALSE

   describe("weak= TRUE black-listing", {
      describe("selectable options - default", {
         it("accepts always good characters", {
            got <- is.safeFileName( alwaysGoodNames, weak= TRUE )
            expect_true( all(got) )
         })
         it("accepts default good characters", {
            got <- is.safeFileName( defaultGoodNames, weak= TRUE )
            expect_true( all(got) )
         })
         it("rejects always bad characters", {
            got <- is.safeFileName( alwaysBadNames, weak= TRUE )
            expect_false( all(got) )
         })
         it("rejects default bad characters", {
            got <- is.safeFileName( defaultBadNames, weak= TRUE )
            expect_false( all(got) )
         })
         it("accepts weakly safe names", {
            got <- is.safeFileName( okWeakName, weak= TRUE )
            expect_true(got)
         })
      }) # END: default selectable
      describe("selectable options - explicit defaults", {
         it("accepts always good characters", {
            got <- is.safeFileName( alwaysGoodNames, weak= TRUE,
               okSpace= TRUE, okTilde= TRUE, okBackslash=FALSE, okComma=FALSE, okColon=FALSE)
            expect_true( all(got) )
         })
         it("accepts default good characters", {
            got <- is.safeFileName( defaultGoodNames, weak= TRUE,
               okSpace= TRUE, okTilde= TRUE, okBackslash=FALSE, okComma=FALSE, okColon=FALSE)
            expect_true( all(got) )
         })
         it("rejects always bad characters", {
            got <- is.safeFileName( alwaysBadNames, weak= TRUE,
               okSpace= TRUE, okTilde= TRUE, okBackslash=FALSE, okComma=FALSE, okColon=FALSE)
            expect_false( all(got) )
         })
         it("rejects default bad characters", {
            got <- is.safeFileName( defaultBadNames, weak= TRUE,
               okSpace= TRUE, okTilde= TRUE, okBackslash=FALSE, okComma=FALSE, okColon=FALSE)
            expect_false( all(got) )
         })
         it("accepts weakly safe names", {
            got <- is.safeFileName( okWeakName, weak= TRUE,
               okSpace= TRUE, okTilde= TRUE, okBackslash=FALSE, okComma=FALSE, okColon=FALSE)
            expect_true(got)
         })
      }) # END: selectable explicit defaults
      describe("selectable options - reverse of defaults", {
         it("Accepts always good characters", {
            got <- is.safeFileName( alwaysGoodNames, weak= TRUE,
               okSpace= FALSE, okTilde= FALSE, okBackslash=TRUE, okComma=TRUE, okColon=TRUE)
            expect_true( all(got) )
         })
         it("accepts default good characters", {
            got <- is.safeFileName( defaultGoodNames, weak= TRUE,
               okSpace= FALSE, okTilde= FALSE, okBackslash=TRUE, okComma=TRUE, okColon=TRUE)
            expect_false( all(got) )
         })
         it("rejects always bad characters", {
            got <- is.safeFileName( alwaysBadNames, weak= TRUE,
               okSpace= FALSE, okTilde= FALSE, okBackslash=TRUE, okComma=TRUE, okColon=TRUE)
            expect_false( all(got) )
         })
         it("rejects default bad characters", {
            got <- is.safeFileName( defaultBadNames, weak= TRUE,
               okSpace= FALSE, okTilde= FALSE, okBackslash=TRUE, okComma=TRUE, okColon=TRUE)
            expect_true( all(got) )
         })
         it("accepts weakly safe names", {
            got <- is.safeFileName( okWeakName, weak= TRUE,
               okSpace= FALSE, okTilde= FALSE, okBackslash=TRUE, okComma=TRUE, okColon=TRUE)
            expect_true(got)
         })

      }) # END: selectable reversed from defaults
      describe("Exists option in context of weak=FALSE",{
         expect_true( file.exists (existingFile ))
         expect_false( file.exists("noSuchFile.exists" ))

         it("doesn't affect default exists option", {
            got <- is.safeFileName( existingFile, weak= TRUE )
            expect_true(got)
            got <- is.safeFileName( "noSuchFile.exists", weak= TRUE)
            expect_true(got)
         })
         it("doesn't affect exists = TRUE", {
            got <- is.safeFileName( existingFile, exists= TRUE, weak= TRUE )
            expect_true(got)
            got <- is.safeFileName( "noSuchFile.exists", exists= TRUE, weak= TRUE)
            expect_false(got)
         })
         it("doesn't affect exists = FALSE", {
            got <- is.safeFileName( existingFile, exists= FALSE, weak= TRUE )
            expect_false(got)
            got <- is.safeFileName( "noSuchFile.exists", exists= FALSE, weak= TRUE)
            expect_true(got)
         })
      }) # END exist=
   }) # END weak=TRUE

   describe("exists= option", {
      expect_true( file.exists (existingFile ))
      expect_false( file.exists("noSuchFile.exists" ))

      it("Both ok by default default exists option", {
         got <- is.safeFileName( existingFile )
         expect_true(got)
         got <- is.safeFileName( "noSuchFile.exists" )
         expect_true(got)
      })
      it("doesn't affect exists = TRUE", {
         got <- is.safeFileName( existingFile, exists= TRUE)
         expect_true(got)
         got <- is.safeFileName( "noSuchFile.exists", exists= TRUE)
         expect_false(got)
      })
      it("doesn't affect exists = FALSE", {
         got <- is.safeFileName( existingFile, exists= FALSE)
         expect_false(got)
         got <- is.safeFileName( "noSuchFile.exists", exists= FALSE)
         expect_true(got)
      })
   }) # END exist= option

}) # END: is.safeFileName()