# SecMML
Machine learning framework supporting secure multiparty computation. Our framework is able to train machine learning models in a privacy preserving manner without leaking information beyond the final model. Any question, please contact **19212010008@fudan.edu.cn**

## Application scenario

There are two practical situations as follow:

1. Several companies each hold their own data sets. They want to train a better model on their union data sets, but would not to leak their own data sets. At first, they share their data to other parties in a secret sharing manner. In this way, each party still has a share of the entire data set. Then, as a party, each company trains the model collaboratively.  When the training is completed, the model is revealed to all parties. Our framework is extensible to support arbitrary number of participants (three+) to train models on the entire data set composed of the data they hold.
2. There are a large number of individual data owners and they do not want their private data to be known by others. Internet companies want to make use of these distributed data to acquire better models. These companies may ﬁrst specify several servers to perform the computation and these servers must be independent of each other. All data owners then send their data to these servers in secret sharing manner. The servers collaboratively train the model with these data and the trained model is ﬁnally revealed to the data owners. The scalability of the framework is that it can support any number of data owners, and any number of servers can be selected as computing parties. 

## BGW protocol

Our framework based on BGW protocol, one of the most fundamental protocol of secure computation presented by Ben-Or,
Goldwasser and Wigderson (BGW) in 1988.

## Machine learning algorithm

Our framework support linear regression, logistic regression, BP neural network and LSTM neural network.

## Key-value models
Non-parametric models are typically based on the key-value operations like get and set. For example, decision tree is constructed node by node and each node represents a specific attribute and corresponding value selection.