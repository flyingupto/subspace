# subspace
这是一个关于论文中描述的 MATLAB 算法的 GitHub 仓库。本仓库包含了实现这个算法的 MATLAB 代码。
## 目录
- [算法简介](#算法简介)
- [开发前的配置要求](#开发前的配置要求)
- [文件目录说明](#文件目录说明)
- [版本控制](#版本控制)
- [作者](#作者)
## 算法简介

这里简要介绍一下论文中描述的算法的背景和主要思想。描述算法的关键概念和原理。

## 开发前的配置要求

MATLAB R2018a及以上

## 安装

1. 克隆或下载本仓库。

2. 打开 MATLAB。

3. 导航到仓库的目录。

## 使用说明

```matlab
%Load a dataset:
[x1,x2,time] = bhs_compare(A, y_noisy,nsamples,burnin, thin);
```

### 文件目录说明

```
filetree 
├── ARCHITECTURE.md
├── LICENSE.txt
├── README.md
│  ├── /variable selection problems/
│  │  ├── bhs_compare.m
│  │  ├── test.m
│  │  └── dataset.mat
└── /util/

```

请阅读**CONTRIBUTING.md** 查阅为该项目做出贡献的开发者。

### 版本控制

该项目使用Git进行版本管理。您可以在repository参看当前可用版本。

### 作者

xxx@xxxx

知乎:xxxx  &ensp; qq:xxxxxx    

 *您也可以在贡献者名单中参看所有参与该项目的开发者。*
