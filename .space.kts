/**
* JetBrains Space Automation
* This Kotlin-script file lets you automate build activities
* For more info, see https://www.jetbrains.com/help/space/automation.html
*/

job("build and test in latest preimage") {
    startOn {
        //push in preimage repo
		gitPush {
            repository = "phantasus-preimage"
            branchFilter {
                +"refs/heads/master"
        	}
    	}
        //push in phantasus repo
        gitPush {
            branchFilter {
            	+"refs/heads/master"
                +"refs/heads/space_jobs"
            }
        }
    }
    container(displayName = "Build R package and copy to share", image = "ctlab.registry.jetbrains.space/p/phantasus/phantasus-containers/phantasus-preimage") {
    	shellScript {
        	content = """
            	R CMD build .
                FILE=${'$'}(ls -1t *.tar.gz | head -n 1)
                cp  ${'$'}FILE $mountDir/share/
            """
        }
    }
    container(displayName = "Check builded package", image = "ctlab.registry.jetbrains.space/p/phantasus/phantasus-containers/phantasus-preimage") {
    	shellScript {
        	content = """
                FILE=${'$'}(ls -1t $mountDir/share/*.tar.gz | head -n 1)
                R CMD check "$mountDir/share/${'$'}FILE"
                bash inst/test_js.sh
            """
        }
    }
    container(displayName = "BioCheck", image = "ctlab.registry.jetbrains.space/p/phantasus/phantasus-containers/phantasus-preimage") {
    	shellScript {
        	content = """
                FILE=${'$'}(ls -1t $mountDir/share/*.tar.gz | head -n 1)
                Rscript -e "library(BiocCheck); BiocCheck(\"${'$'}{FILE}\")"
                Rscript -e 'covr::codecov()'
            """
        }
    }
    container(displayName = "code covr", image = "ctlab.registry.jetbrains.space/p/phantasus/phantasus-containers/phantasus-preimage") {
    	shellScript {
        	content = """
                Rscript -e "install.packages(\"covr\")"
                Rscript -e 'covr::codecov()'
            """
        }
    }
    
}

