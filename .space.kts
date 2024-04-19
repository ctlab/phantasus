/**
* JetBrains Space Automation
* This Kotlin-script file lets you automate build activities
* For more info, see https://www.jetbrains.com/help/space/automation.html
*/

job("build and test in latest preimage") {
    startOn {
        //push in phantasus repo
        gitPush {
           anyBranchMatching {
                +"master"
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
                        export _R_CHECK_FORCE_SUGGESTS_=FALSE
                        R CMD check "${'$'}FILE" --no-manual
                    """
                }
            }
        }
        container(displayName = "Check js", image = "ctlab.registry.jetbrains.space/p/phantasus/phantasus-containers/phantasus-preimage") {
            shellScript {
                content = """
		    add-apt-repository -y  ppa:mozillateam/ppa
		    echo '
Package: *
Pin: release o=LP-PPA-mozillateam
Pin-Priority: 1001
' | tee /etc/apt/preferences.d/mozilla-firefox
                    apt-get update && apt-get install -y --no-install-recommends nodejs npm firefox
                    bash inst/test_js.sh
                """
            }
    	}
    }


    kaniko {
    	beforeBuildScript {
            // Create an env variable BRANCH,
            // use env var to get full branch name,
            // leave only the branch name without the 'refs/heads/' path
            content = """
                export BRANCH=${'$'}(echo ${'$'}JB_SPACE_GIT_BRANCH | cut -d'/' -f 3)
                export DOCKER_TAG=${'$'}BRANCH-${'$'}JB_SPACE_EXECUTION_NUMBER
                case ${'$'}BRANCH in 
                  "master") export LATEST="latest" ;;
                  *)  export LATEST="test" ;;
                esac
            """
        }  
        build {
        	context = "."
            dockerfile = "Dockerfile"
            args["PREIMAGE_NAME"] = "ctlab.registry.jetbrains.space/p/phantasus/phantasus-containers/phantasus-preimage"
            args["TARGET_BRANCH"] = "\$BRANCH"
            args["PHANTASUS_BUILD"] = "\$DOCKER_TAG"
            labels["vendor"] = "ctlab"
        }
         push("ctlab.registry.jetbrains.space/p/phantasus/phantasus-containers/phantasus") {
            tags("\$DOCKER_TAG", "\$LATEST")
        }
    }
    
}

