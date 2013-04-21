################################################################################
## topic-networks.R A sample analysis of a document collection using some     ##
##    standard NLP techniques and a new approach with topic networks.         ##
##                                                                            ##
## To run this code you must first download and extract the sample corpus     ##
##  corpus.zip                                                                ##
## and set the corpus.directory variable to point to the corpus directory.    ##
##                                                                            ##
## Requires the following packages:                                           ##
##    igraph                                                                  ##
##    ggplot2                                                                 ##
##    RColorBrewer                                                            ##
##    slam                                                                    ##
##    tm                                                                      ##
##    topicmodels                                                             ##
##    wordcloud                                                               ##
##                                                                            ##
################################################################################
# Set the path to the corpus directory.                                        #
# ---------------------------------------------------------------------------- #
corpus.directory <- '/home/birksworks/Projects/topic-networks/corpus'
# ---------------------------------------------------------------------------- #
# Load required libraries.                                                     #
# ---------------------------------------------------------------------------- #
require(igraph, quietly=TRUE)
require(RColorBrewer, quietly=TRUE)
require(ggplot2, quietly=TRUE)
require(slam, quietly=TRUE)
require(tm, quietly=TRUE)
require(topicmodels, quietly=TRUE)
require(wordcloud, quietly=TRUE)
# ---------------------------------------------------------------------------- #
# Load the corpus data.                                                        #
# ---------------------------------------------------------------------------- #
corpus.source <- DirSource(corpus.directory, encoding="UTF-8")
corpus <- Corpus(corpus.source)
corpus.copy <-corpus
# ---------------------------------------------------------------------------- #
# Build the document term matrix.                                              #
# ---------------------------------------------------------------------------- #
corpus.dtm <- DocumentTermMatrix(
  corpus, 
  control = list(
    stemming = TRUE, 
    stopwords = TRUE, 
    minWordLength = 3,
    removeNumbers = TRUE, 
    removePunctuation = TRUE
  )
)
# ---------------------------------------------------------------------------- #
# Visualize frequent words.                                                    #
# ---------------------------------------------------------------------------- #
word.freq <- sort(col_sums(corpus.dtm), decreasing=TRUE)
# barplot with the top 20 most frequent terms
top.terms <- head(word.freq, n=20)
# complete the stems and fix missing values and errors:
completions <- stemCompletion(names(top.terms), dictionary=corpus.copy, type="prevalent")
completions<-ifelse(completions == "", names(top.terms), completions)
names(top.terms) <-completions
names(top.terms)[10] <-"every"
pdf("mostfreq.pdf")
op <- par(mar = c(4,6.1,.1,.2))
barplot(top.terms, las=2, horiz=TRUE)
dev.off()
# wordcloud
pal2 <- brewer.pal(8,"Dark2")
top.terms <- head(word.freq, n=1200)
completions <- stemCompletion(names(top.terms), dictionary=corpus.copy, type="prevalent")
completions<-ifelse(completions == "", names(top.terms), completions)
names(top.terms) <-completions
pdf("wordcloud.pdf")
par(mar = c(0,0,0,0))
wordcloud(words=names(top.terms), freq=top.terms, min.freq=20, random.order=F, colors=pal2, rot.per=.15)
dev.off()
# ---------------------------------------------------------------------------- #
# Build the tf-idf matrix and use it to filter the document term matrix.       #
# ---------------------------------------------------------------------------- #
term.tfidf <- tapply(corpus.dtm$v/row_sums(corpus.dtm)[corpus.dtm$i], corpus.dtm$j, mean) * log2(nDocs(corpus.dtm)/col_sums(corpus.dtm > 0))
dtm <- corpus.dtm[,term.tfidf >= 0.01]
dtm <- dtm[row_sums(dtm) > 0,]
# ---------------------------------------------------------------------------- #
# Build a topic model and collect relevant data in data frames.                #
# ---------------------------------------------------------------------------- #
corpus.tm <- LDA(dtm, k = 22)
corpus.tm.terms <- terms(corpus.tm, 3)
corpus.tm.topics <- topics(corpus.tm, 4)
topic.labels <- apply(corpus.tm.terms, 2, function(x) paste(x, collapse=", "))
document.labels <- colnames(corpus.tm.topics)
dt.df <- data.frame(
  document=rep(document.labels, each=4),
  topic=as.vector(corpus.tm.topics)
)
dt.df$document <- as.numeric(gsub(".txt","",dt.df$document))
dt.df <- dt.df[order(dt.df$document),]
# ---------------------------------------------------------------------------- #
# Build the incidence matrice and document-topic network.                      #
# ---------------------------------------------------------------------------- #
dt.matrix <- as.matrix(table(dt.df))
dt.network <- graph.incidence(dt.matrix)
# ---------------------------------------------------------------------------- #
# Visualize the document-topic network.                                        #
# ---------------------------------------------------------------------------- #
n.docs <- nrow(dt.matrix)
n.topics <- ncol(dt.matrix)
n.vertices <- n.docs + n.topics
V(dt.network)$color[1:n.docs] <- rgb(1,0,0,.4)
V(dt.network)$color[(n.docs+1):n.vertices] <- rgb(0,1,0,.5)
V(dt.network)$label <- NA
V(dt.network)$size[1:n.docs] <- 2
V(dt.network)$size[(n.docs+1):n.vertices] <- 6
E(dt.network)$width <- .5
E(dt.network)$color <- rgb(.5,.5,0,.4)
png("dtn.png")
plot(dt.network, layout=layout.fruchterman.reingold)
dev.off()
# ---------------------------------------------------------------------------- #
# Create the topic network.                                                    #
# ---------------------------------------------------------------------------- #
topic.network.matrix <- t(dt.matrix) %*% dt.matrix
topic.network <- graph.adjacency(topic.network.matrix, mode = "undirected")
E(topic.network)$weight <- count.multiple(topic.network)
topic.network <- simplify(topic.network)
# Set vertex attributes
V(topic.network)$label <- topic.labels
V(topic.network)$label.color <- rgb(0,0,0,1)
V(topic.network)$label.cex <- .75
V(topic.network)$size <- 8
#V(topic.network)$frame.color <- NA
V(topic.network)$color <- rgb(0,1,0,.6)

# Set edge gamma according to edge weight
egam <- (log(E(topic.network)$weight)+.3)/max(log(E(topic.network)$weight)+.3)
E(topic.network)$color <- rgb(.5,.5,0,egam)
pdf("topic-network.pdf")
plot(topic.network, layout=layout.kamada.kawai)
dev.off()
# ---------------------------------------------------------------------------- #
# Topic overlap.                                                               #
# ---------------------------------------------------------------------------- #
overlap.matrix <- topic.network.matrix / diag(topic.network.matrix)
overlap.network <- graph.adjacency(overlap.matrix, weighted=T)
# Degree
V(overlap.network)$degree <- degree(overlap.network)
# Betweenness centrality
V(overlap.network)$btwcnt <- betweenness(overlap.network)
# Plot the connection strength:
plot(topic.network, layout=layout.kamada.kawai)
pdf("density.pdf")
plot(density(overlap.matrix), main=NA, xlab=NA)
dev.off()
overlap.matrix[overlap.matrix < 0.2] <- 0
overlap.network <- graph.adjacency(overlap.matrix, weighted=T)
overlap.network <- simplify(overlap.network, remove.multiple=FALSE, remove.loops=TRUE)
overlap.network$layout <- layout.kamada.kawai(overlap.network)
V(overlap.network)$label <- topic.labels
tkplot(overlap.network)
overlap.network$layout <- tkplot.getcoords(1)

# Set vertex attributes
V(overlap.network)$label.color <- rgb(0,0,.2,.6)
V(overlap.network)$label.cex <- .75
V(overlap.network)$size <- 6
#V(topic.network)$frame.color <- NA
V(overlap.network)$color <- rgb(0,1,0,.6)

# Set edge gamma according to edge weight
egam <- (E(overlap.network)$weight+.1)/max(E(overlap.network)$weight+.1)
E(overlap.network)$color <- rgb(.5,.5,0,egam)
E(overlap.network)$arrow.size <- .3
V(overlap.network)$label.cex <- degree(overlap.network)/(max(degree(overlap.network)/2))+ .3
pdf("overlap-network.pdf")
plot(overlap.network)
dev.off()
