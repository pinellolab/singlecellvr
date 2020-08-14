AFRAME.registerComponent('loading', {
    schema: {
        time: {type: 'number', default: 5},
        show: {type: 'boolean', default: true}
    },
    init: function() {
        this.initialized = performance.now();
        var geometry = new THREE.PlaneBufferGeometry( .5, .5 );
        var material = new THREE.MeshBasicMaterial( {color: 0xffff00, side: THREE.DoubleSide} );
        var plane = new THREE.Mesh( geometry, material );
        this.el.object3D.add(plane);
    },
    update: function(oldData) {
        if (this.data.show !== oldData.show && !this.data.show) {
            const dismissTime = performance.now();
            const elapsedSeconds = (dismissTime - this.initialized) / 1000;
            if (elapsedSeconds > this.data.time) {
                this.el.object3D.visible = false;
            } else {
                setTimeout(() => {
                    this.el.object3D.visible = false;
                }, ((this.data.time - elapsedSeconds) * 1000))
            }
        }
    }

});