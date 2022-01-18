/**
* JetBrains Space Automation
* This Kotlin-script file lets you automate build activities
* For more info, see https://www.jetbrains.com/help/space/automation.html
*/

job("build and test in latest preimage") {
    startOn {
        //push in phantasus repo
        gitPush {
           branchFilter {
                +"refs/heads/master"
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
    parallel{
    	sequential{
            container(displayName = "Check builded package", image = "ctlab.registry.jetbrains.space/p/phantasus/phantasus-containers/phantasus-preimage") {
                shellScript {
                    content = """
                			FILE=${'$'}(ls -1t $mountDir/share/*.tar.gz | head -n 1)
                			R CMD check "${'$'}FILE" --no-manual
                    """
                }
            }
            container(displayName = "BioCheck", image = "ctlab.registry.jetbrains.space/p/phantasus/phantasus-containers/phantasus-preimage") {
                shellScript {
                    content = """
                			FILE=${'$'}(ls -1t $mountDir/share/*.tar.gz | head -n 1)
                            Rscript -e "BiocManager::install(\"BiocCheck\")"
                            Rscript -e "library(BiocCheck); BiocCheck(\"${'$'}{FILE}\")"
                    """
                }
            }
            container(displayName = "code covr", image = "ctlab.registry.jetbrains.space/p/phantasus/phantasus-containers/phantasus-preimage") {
                shellScript {
                    content = """
                            Rscript -e "install.packages(\"covr\")"
                            Rscript -e 'covr::package_coverage(quiet = FALSE)'
                     """
                }
            }
        }
        container(displayName = "Check js", image = "ctlab.registry.jetbrains.space/p/phantasus/phantasus-containers/phantasus-preimage") {
            shellScript {
                content = """
                    apt-get update && apt-get install -y --no-install-recommends nodejs npm firefox
                    bash inst/test_js.sh
                """
            }
    	}
    }

    
    

    docker {
    	beforeBuildScript {
            // Create an env variable BRANCH,
            // use env var to get full branch name,
            // leave only the branch name without the 'refs/heads/' path
            content = """
                export BRANCH=${'$'}(echo ${'$'}JB_SPACE_GIT_BRANCH | cut -d'/' -f 3)
                export DOCKER_TAG=${'$'}BRANCH-${'$'}JB_SPACE_EXECUTION_NUMBER
                export GITHUB_PAT=${'$'}JB_SPACE_API_URL/${'$'}JB_SPACE_PROJECT_KEY
                case ${'$'}BRANCH in 
                  "master") export LATEST="latest" ;;
                  *)  export LATEST="test" ;;
                esac
            """
        }  
        build {
        	context = "."
            file = "Dockerfile"
            args["PREIMAGE_NAME"] = "ctlab.registry.jetbrains.space/p/phantasus/phantasus-containers/phantasus-preimage"
            args["TARGET_BRANCH"] = "\$BRANCH"
            args["PHANTASUS_BUILD"] = "\$DOCKER_TAG"
            args["GITHUB_PAT"] = "\$GITHUB_PAT"
            labels["vendor"] = "ctlab"
        }
         push("ctlab.registry.jetbrains.space/p/phantasus/phantasus-containers/phantasus") {
            tags("\$DOCKER_TAG", "\$LATEST")
        }
    }
    
}

