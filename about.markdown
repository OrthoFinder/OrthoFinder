---
layout: page
title: About
permalink: /about/
---


## About

### Why OrthoFinder?

In the coming decades, millions of species will have their genomes sequenced. This data has the power to transform comparative genomics, but we need scalable tools to handle this volume of data. OrthoFinder meets this challenge, with a new scalable framework that allows rapid orthology inference on thousands of species, without compromising accuracy. 
We focus on:

<ul class="tight-bullets">
  <li>
    <div class="bullet-header">
      <span class="bullet">•</span>
      <span class="bullet-title">Accuracy</span>
    </div>
    <div class="bullet-text">
      OrthoFinder is the overall best performing tool across orthogroup and ortholog benchmarking.
    </div>
  </li>
  <li>
    <div class="bullet-header">
      <span class="bullet">•</span>
      <span class="bullet-title">Scalability</span>
    </div>
    <div class="bullet-text">
      OrthoFinder allows you to analyse thousands of species on conventional computing resources.
    </div>
  </li>
  <li>
    <div class="bullet-header">
      <span class="bullet">•</span>
      <span class="bullet-title">Usability</span>
    </div>
    <div class="bullet-text">
      OrthoFinder only requires a single fasta file for each species, and outputs a vast array of comparative genomic information including gene trees and gene duplication events.
    </div>
  </li>
  <li>
    <div class="bullet-header">
      <span class="bullet">•</span>
      <span class="bullet-title">Flexibility</span>
    </div>
    <div class="bullet-text">
      Users can easily swap to different tree-inference and multiple sequence alignment tools, thanks to OrthoFinder’s modular design.
    </div>
  </li>
</ul>

### Development and community
OrthoFinder is developed by [Steve Kelly’s research group](http://www.stevekellylab.com/) at the University of Oxford. We have an active community of users and developers on github, who can provide help and assistance 
 
### Meet&nbsp;the&nbsp;Team {#team}

#### Current&nbsp;Team&nbsp;Member {#team-current}

<div class="team-grid">
{% assign current_team = site.data.team | where: "status", "Current" %}
{% for member in current_team %}
  {% include team-card.html member=member %}
{% endfor %}
</div>

---

#### Former&nbsp;Team&nbsp;Member {#team-former}

<div class="team-grid">
{% assign former_team = site.data.team | where: "status", "Former" %}
{% for member in former_team %}
  {% include team-card.html member=member %}
{% endfor %}
</div>
