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
    container(displayName = "", image = "ctlab.registry.jetbrains.space/p/phantasus/phantasus-containers/phantasus-preimage:0.2-using_ocpu") {
    	shellScript {
        	content = """
            	R CMD build .
                FILE=$(ls -1t *.tar.gz | head -n 1)
                R CMD check "${'$'}FILE"
                bash inst/test_js.sh
                Rscript -e "library(BiocCheck); BiocCheck(\"${'$'}{FILE}\")"
                Rscript -e 'covr::codecov()'
            """
        }
    }
}

