pipeline
{
    agent any 
    
    stages
    {
        stage ("identification")
        {
            steps
            {
                echo "Committer: ${env.GIT_COMMITTER_NAME}."
                echo "Email:     ${env.GIT_COMMITTER_EMAIL}."
                echo "Branch:    ${env.BRANCH_NAME}."
            }
        }

        stage ("build")
        {
            steps
            {
                echo "Building . . ."
                sh 'make pi'
                sh 'make e'
            }
        }

        stage ("deploy")
        {
            steps
            {
                echo "Deploying . . ."
                fileOperations([fileCopyOperation(
                    excludes: '',
                    flattenFiles: false,
                    includes: '*.exe',
                    targetLocation: "/home/nmoulton/Documents/JenkinsBuilds"
                )])

                //echo "Contents of output folder:"
                //sh 'ls /home/nmoulton/Documents/JenkinsBuilds/'
            }
        }
    }
}